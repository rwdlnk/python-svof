# cg_pressure.py — Jacobi-preconditioned CG pressure solver for NMAT=2
#
# Replaces the inherently sequential Gauss-Seidel PRESSIT with a fully
# parallel conjugate gradient solver. Matrix-free: the Laplacian (SpMV)
# is applied via a Numba kernel using the variable-coefficient 5-point stencil.
#
# The pressure Poisson equation solved:
#
#   sum_faces [ 2*dt / rho_face * (dp_neighbor - dp_ij) / dx_face ] = div(u_tilde)
#
# where div(u_tilde) = rdx*(U[i,j]-U[i-1,j]) + rdy*(V[i,j]-V[i,j-1])
#
# This is A*dp = b, where A is the 5-point variable-coefficient Laplacian
# and b = div(u_tilde). After solving, we update:
#   P += dp
#   U[i,j] += 2*dt*comg*dp / rhoxr  (right face)
#   U[i-1,j] -= 2*dt*comg*dp / rhoxl (left face)
#   etc. for V
import numpy as np
from numba import njit


EMF = 1.0e-6

# =====================================================================
# Numba kernels
# =====================================================================

@njit(cache=True)
def _compute_rhs(U, V, F, BETA, rdx, rdy, rxi,
                 i_lo, i_hi, j_lo, j_hi, cyl, nmat, rhs):
    """
    Compute RHS = div(u_tilde) for the pressure Poisson equation.
    """
    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                rhs[i, j] = 0.0
                continue
            if nmat != 2 and F[i, j] < EMF:
                rhs[i, j] = 0.0
                continue

            dij = rdx[i] * (U[i, j] - U[i - 1, j]) + rdy[j] * (V[i, j] - V[i, j - 1])
            if cyl > 0.5:
                nrxi = rxi.shape[0]
                if i < nrxi:
                    dij += cyl * rxi[i] * (U[i, j] + U[i - 1, j]) / 2.0

            rhs[i, j] = dij


@njit(cache=True)
def _compute_coeffs(F, BETA, delx, dely, rdx, rdy,
                    i_lo, i_hi, j_lo, j_hi,
                    dt, rhofc, rhod, comg, rcsq,
                    ar_arr, al_arr, at_arr, ab_arr, diag_arr):
    """
    Precompute the stencil coefficients for the pressure Laplacian.

    The coefficient for (e.g.) the right neighbor is:
      ar = 2*dt*comg * rdx[i] / rhoxr

    where rhoxr = (rho_ij * dx_ip + rho_ip * dx_i)

    These match exactly what pressit uses for DPTC velocity correction.
    The diagonal is minus the sum of off-diagonals (plus compressibility).
    """
    ndelx = delx.shape[0]
    ndely = dely.shape[0]

    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                ar_arr[i, j] = 0.0
                al_arr[i, j] = 0.0
                at_arr[i, j] = 0.0
                ab_arr[i, j] = 0.0
                diag_arr[i, j] = 1.0  # identity for obstacles
                continue

            dx_i = delx[i]
            dx_ip = delx[min(i + 1, ndelx - 1)]
            dx_im = delx[max(i - 1, 0)]
            dy_j = dely[j]
            dy_jp = dely[min(j + 1, ndely - 1)]
            dy_jm = dely[max(j - 1, 0)]

            rho_ij = rhofc + rhod * F[i, j]

            # Right face
            rhoxr = rho_ij * dx_ip + (rhofc + rhod * F[i + 1, j]) * dx_i
            ar = 2.0 * dt * comg * rdx[i] / rhoxr if rhoxr != 0.0 else 0.0

            # Left face
            rhoxl = (rhofc + rhod * F[i - 1, j]) * dx_i + rho_ij * dx_im
            al = 2.0 * dt * comg * rdx[i] / rhoxl if rhoxl != 0.0 else 0.0

            # Top face
            rhoyt = rho_ij * dy_jp + (rhofc + rhod * F[i, j + 1]) * dy_j
            at = 2.0 * dt * comg * rdy[j] / rhoyt if rhoyt != 0.0 else 0.0

            # Bottom face
            rhoyb = (rhofc + rhod * F[i, j - 1]) * dy_j + rho_ij * dy_jm
            ab = 2.0 * dt * comg * rdy[j] / rhoyb if rhoyb != 0.0 else 0.0

            # Block neighbors that are obstacles
            if BETA[i + 1, j] < 0.0:
                ar = 0.0
            if BETA[i - 1, j] < 0.0:
                al = 0.0
            if BETA[i, j + 1] < 0.0:
                at = 0.0
            if BETA[i, j - 1] < 0.0:
                ab = 0.0

            ar_arr[i, j] = ar
            al_arr[i, j] = al
            at_arr[i, j] = at
            ab_arr[i, j] = ab

            # Diagonal = -(ar + al + at + ab) - compressibility term
            comp = 0.0
            if rcsq != 0.0:
                rhor = (rhofc + rhod) / rho_ij
                comp = rhor * rcsq * comg / dt
            diag_arr[i, j] = -(ar + al + at + ab) - comp


@njit(cache=True)
def _apply_laplacian(p, ar, al, at, ab, diag, BETA,
                     i_lo, i_hi, j_lo, j_hi, Ap):
    """
    Matrix-free SpMV: Ap = A * p using precomputed coefficients.
    A*p_ij = diag*p_ij + ar*p_{i+1,j} + al*p_{i-1,j} + at*p_{i,j+1} + ab*p_{i,j-1}
    """
    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                Ap[i, j] = p[i, j]
                continue

            Ap[i, j] = (diag[i, j] * p[i, j] +
                        ar[i, j] * p[i + 1, j] +
                        al[i, j] * p[i - 1, j] +
                        at[i, j] * p[i, j + 1] +
                        ab[i, j] * p[i, j - 1])


@njit(cache=True)
def _local_dot(a, b, BETA, i_lo, i_hi, j_lo, j_hi):
    """Local dot product over interior fluid cells."""
    s = 0.0
    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                continue
            s += a[i, j] * b[i, j]
    return s


@njit(cache=True)
def _local_max_abs(a, BETA, i_lo, i_hi, j_lo, j_hi):
    """Local infinity norm over interior fluid cells."""
    mx = 0.0
    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                continue
            mx = max(mx, abs(a[i, j]))
    return mx


@njit(cache=True)
def _axpy(y, alpha, x, BETA, i_lo, i_hi, j_lo, j_hi):
    """y += alpha * x over interior fluid cells."""
    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                continue
            y[i, j] += alpha * x[i, j]


@njit(cache=True)
def _precond_apply(z, diag, r, BETA, i_lo, i_hi, j_lo, j_hi):
    """z = (1/diag) * r — Jacobi preconditioner."""
    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                z[i, j] = 0.0
                continue
            d = diag[i, j]
            z[i, j] = r[i, j] / d if abs(d) > 1.0e-30 else 0.0


@njit(cache=True)
def _update_search(p_cg, z, beta_cg, BETA, i_lo, i_hi, j_lo, j_hi):
    """p = z + beta * p (update search direction)."""
    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                continue
            p_cg[i, j] = z[i, j] + beta_cg * p_cg[i, j]


@njit(cache=True)
def _apply_pressure_correction(U, V, P, dp, F, BETA,
                               delx, dely,
                               i_lo, i_hi, j_lo, j_hi,
                               dt, rhofc, rhod, comg):
    """
    After CG solves for dp, update P += dp and correct U, V.

    Each cell corrects only its RIGHT face U[i,j] and TOP face V[i,j],
    applying the combined correction from both cells sharing that face:
      U[i,j] += 2*dt*comg * (dp[i,j] - dp[i+1,j]) / rho_face

    This avoids writing to ghost-cell velocities, which is critical for
    MPI: writing to a ghost and then halo_exchanging would overwrite the
    correction, losing the neighbor cell's contribution.

    Requires dp ghost cells to be filled (halo_exchange + Neumann) before
    calling this function.
    """
    ndelx = delx.shape[0]
    ndely = dely.shape[0]

    for i in range(i_lo, i_hi):
        for j in range(j_lo, j_hi):
            if BETA[i, j] < 0.0:
                continue

            P[i, j] += dp[i, j]

            coeff = 2.0 * dt * comg
            dx_i = delx[i]
            dx_ip = delx[min(i + 1, ndelx - 1)]
            dy_j = dely[j]
            dy_jp = dely[min(j + 1, ndely - 1)]

            # Right face U[i,j]: combined correction from cells i and i+1
            # rho_face = rho_i * dx_{i+1} + rho_{i+1} * dx_i
            # (same as rhoxr for cell i and rhoxl for cell i+1)
            if BETA[i + 1, j] >= 0.0:
                rho_face = ((rhofc + rhod * F[i, j]) * dx_ip +
                            (rhofc + rhod * F[i + 1, j]) * dx_i)
                if rho_face != 0.0:
                    U[i, j] += coeff * (dp[i, j] - dp[i + 1, j]) / rho_face

            # Top face V[i,j]: combined correction from cells j and j+1
            if BETA[i, j + 1] >= 0.0:
                rho_face = ((rhofc + rhod * F[i, j]) * dy_jp +
                            (rhofc + rhod * F[i, j + 1]) * dy_j)
                if rho_face != 0.0:
                    V[i, j] += coeff * (dp[i, j] - dp[i, j + 1]) / rho_face


# =====================================================================
# Helper: Neumann ghost BCs for CG arrays
# =====================================================================

def _apply_neumann_ghosts(arr, decomp):
    """Copy interior boundary values to physical-wall ghost cells
    so that dp/dn = 0 (Neumann BC), matching GS solver behavior."""
    if decomp.has_left_wall:
        arr[0, :] = arr[1, :]
    if decomp.has_right_wall:
        arr[-1, :] = arr[-2, :]
    if decomp.has_bottom_wall:
        arr[:, 0] = arr[:, 1]
    if decomp.has_top_wall:
        arr[:, -1] = arr[:, -2]


# =====================================================================
# Python CG wrapper with MPI
# =====================================================================

def cg_pressure_solve(decomp, mesh, fields, state, bc_func,
                      max_iter=500, verbose=False):
    """
    Jacobi-preconditioned Conjugate Gradient pressure solver.

    Solves: A * dp = b where
      b = div(u_tilde)
      A_ij = -(ar+al+at+ab) * dp_ij + ar*dp_{i+1} + al*dp_{i-1} + ...
    with ar = 2*dt*comg*rdx/rhoxr matching SOLA's velocity correction.

    Then applies P += dp and velocity corrections.
    """
    from mpi4py import MPI
    from .halo import halo_exchange

    comm = decomp.comm
    im1 = mesh.im1
    jm1 = mesh.jm1
    i_lo, i_hi = 1, im1
    j_lo, j_hi = 1, jm1
    shape = (mesh.imax, mesh.jmax)

    dt = state.delt
    epsi = state.epsi
    rhofc = state.rhofc
    rhod = state.rhof - state.rhofc
    csq = state.csq
    rcsq = 1.0 / (state.rhof * csq) if csq > 0.0 else 0.0
    cyl = state.cyl

    ctos = dt * state.rdtexp
    comg = min(ctos * ctos, 1.0)

    # Allocate working arrays
    rhs = np.zeros(shape, dtype=float)
    Ap = np.zeros(shape, dtype=float)
    r = np.zeros(shape, dtype=float)
    z = np.zeros(shape, dtype=float)
    p_cg = np.zeros(shape, dtype=float)
    dp = np.zeros(shape, dtype=float)

    # Precompute stencil coefficients
    ar = np.zeros(shape, dtype=float)
    al = np.zeros(shape, dtype=float)
    at = np.zeros(shape, dtype=float)
    ab = np.zeros(shape, dtype=float)
    diag = np.zeros(shape, dtype=float)

    _compute_coeffs(fields.F, fields.BETA,
                    mesh.delx, mesh.dely, mesh.rdx, mesh.rdy,
                    i_lo, i_hi, j_lo, j_hi,
                    dt, rhofc, rhod, comg, rcsq,
                    ar, al, at, ab, diag)

    # Compute RHS = div(u_tilde)
    _compute_rhs(fields.U, fields.V, fields.F, fields.BETA,
                 mesh.rdx, mesh.rdy, mesh.rxi,
                 i_lo, i_hi, j_lo, j_hi,
                 cyl, state.nmat, rhs)

    # Initial residual: r = b - A*dp (dp=0 initially, so r = b = rhs)
    r[:] = rhs

    # Apply preconditioner: z = M^{-1} * r
    _precond_apply(z, diag, r, fields.BETA, i_lo, i_hi, j_lo, j_hi)

    # Initial search direction
    p_cg[:] = z

    # rz_old = <r, z>
    rz_local = _local_dot(r, z, fields.BETA, i_lo, i_hi, j_lo, j_hi)
    rz_old = comm.allreduce(rz_local, op=MPI.SUM)

    state.flg = 1.0
    state.iter = 0
    res_max = 1.0e30

    for iteration in range(max_iter):
        state.iter = iteration + 1

        # Halo exchange on search direction
        halo_exchange(decomp, p_cg)

        # Neumann BC on search direction at physical walls:
        # p_cg[ghost] = p_cg[interior] gives dp/dn = 0, matching the
        # GS solver's behavior where BC sets P[ghost] = P[interior]
        # each sweep.  Without this, ghost = 0 gives Dirichlet BC
        # which causes wrong pressure corrections near walls.
        _apply_neumann_ghosts(p_cg, decomp)

        # Matrix-vector product: Ap = A * p_cg
        _apply_laplacian(p_cg, ar, al, at, ab, diag, fields.BETA,
                         i_lo, i_hi, j_lo, j_hi, Ap)

        # alpha = rz_old / <p, Ap>
        pAp_local = _local_dot(p_cg, Ap, fields.BETA, i_lo, i_hi, j_lo, j_hi)
        pAp = comm.allreduce(pAp_local, op=MPI.SUM)

        if abs(pAp) < 1.0e-30:
            break

        alpha_cg = rz_old / pAp

        # dp += alpha * p_cg
        _axpy(dp, alpha_cg, p_cg, fields.BETA, i_lo, i_hi, j_lo, j_hi)

        # r -= alpha * Ap
        _axpy(r, -alpha_cg, Ap, fields.BETA, i_lo, i_hi, j_lo, j_hi)

        # Check convergence: max|r| (this is max|divergence| after correction)
        res_local = _local_max_abs(r, fields.BETA, i_lo, i_hi, j_lo, j_hi)
        res_max = comm.allreduce(res_local, op=MPI.MAX)

        if res_max < epsi:
            state.flg = 0.0
            break

        # Apply preconditioner: z = M^{-1} * r
        _precond_apply(z, diag, r, fields.BETA, i_lo, i_hi, j_lo, j_hi)

        # beta = <r_new, z_new> / <r_old, z_old>
        rz_new_local = _local_dot(r, z, fields.BETA, i_lo, i_hi, j_lo, j_hi)
        rz_new = comm.allreduce(rz_new_local, op=MPI.SUM)

        if abs(rz_old) < 1.0e-30:
            break

        beta_cg = rz_new / rz_old
        rz_old = rz_new

        # p = z + beta * p
        _update_search(p_cg, z, beta_cg, fields.BETA, i_lo, i_hi, j_lo, j_hi)

    else:
        state.fnoc = 1.0
        if decomp.rank == 0:
            print(f"CG: no convergence after {max_iter} iterations "
                  f"(res={res_max:.3e})")

    # Exchange dp so each rank has neighbor dp values in ghost cells.
    # This lets _apply_pressure_correction compute the combined correction
    # for boundary faces: U[i,j] += coeff*(dp[i,j] - dp[i+1,j])/rho_face
    # where dp[i+1,j] comes from the neighboring rank via halo.
    halo_exchange(decomp, dp)
    _apply_neumann_ghosts(dp, decomp)

    # Apply pressure correction to P, U, V
    _apply_pressure_correction(fields.U, fields.V, fields.P, dp,
                               fields.F, fields.BETA,
                               mesh.delx, mesh.dely,
                               i_lo, i_hi, j_lo, j_hi,
                               dt, rhofc, rhod, comg)

    # Halo exchange on corrected fields
    halo_exchange(decomp, fields.P, fields.U, fields.V)

    # Apply physical BCs
    if bc_func is not None:
        bc_func()

    # Diagnostic: check actual divergence after correction + BCs
    if verbose:
        div_check = np.zeros(shape, dtype=float)
        _compute_rhs(fields.U, fields.V, fields.F, fields.BETA,
                      mesh.rdx, mesh.rdy, mesh.rxi,
                      i_lo, i_hi, j_lo, j_hi,
                      cyl, state.nmat, div_check)
        local_div_max = _local_max_abs(div_check, fields.BETA, i_lo, i_hi, j_lo, j_hi)
        global_div_max = comm.allreduce(local_div_max, op=MPI.MAX)
        if decomp.rank == 0:
            print(f"CG: {state.iter} iters, res={res_max:.3e}, "
                  f"post-correction div_max={global_div_max:.3e}")
