# solver.py
import numpy as np
from numba import njit
from .mesh import Mesh
from .fields import Fields
from .state import RunState


def tilde_step(mesh: Mesh, fields: Fields, state: RunState) -> None:
    """
    Vectorized Python version of Fortran SUBROUTINE TILDE for U and V.
    Computes provisional velocities from advection, pressure gradient,
    viscosity, and gravity using the WUDS (alpha) scheme.
    """

    UN, VN = fields.UN, fields.VN
    P = fields.P
    F = fields.F

    delx = mesh.delx
    dely = mesh.dely
    rdx = mesh.rdx
    rdy = mesh.rdy
    rx = mesh.rx

    dt = state.delt
    gx, gy = state.gx, state.gy
    vnu = state.vnu
    vnuc = state.vnuc

    rhof = state.rhof
    rhofc = state.rhofc
    rhod = rhof - rhofc
    alpha = state.alpha

    im1, jm1 = mesh.im1, mesh.jm1

    # Interior range matching old scalar code: i=1..im1-2, j=1..jm1-2
    # This is the safe range where all stencil accesses (i-1, i+1, j-1, j+1)
    # and mesh arrays (delx[i+1], dely[j+1]) stay in bounds.
    si = slice(1, im1 - 1)   # i indices
    sj = slice(1, jm1 - 1)   # j indices

    # ======== U momentum ========
    # Interior i=1..im1-2: ni = im1-2 elements
    # Interior j=1..jm1-2: nj = jm1-2 elements
    # delx[i+1] when i=1..im1-2 means delx[2..im1-1]
    dx_i  = delx[si, np.newaxis]              # delx[i]   (ni, 1)
    dx_ip = delx[2:im1, np.newaxis]           # delx[i+1]
    rx_i  = rdx[si, np.newaxis]               # rdx[i]
    rx_ip = rdx[2:im1, np.newaxis]            # rdx[i+1]

    dy_j  = dely[np.newaxis, sj]              # dely[j]   (1, nj)
    dy_jp = dely[np.newaxis, 2:jm1]           # dely[j+1]
    dy_jm = dely[np.newaxis, 0:jm1-2]        # dely[j-1]
    ry_j  = rdy[np.newaxis, sj]               # rdy[j]
    ry_jp = rdy[np.newaxis, 2:jm1]            # rdy[j+1]
    ry_jm = rdy[np.newaxis, 0:jm1-2]         # rdy[j-1]

    rdelx = 1.0 / (dx_i + dx_ip)

    # Extract UN stencil values for i=1..im1-2, j=1..jm1-2
    ni = im1 - 2   # number of interior i points
    nj = jm1 - 2   # number of interior j points
    un_ij   = UN[si, sj]
    un_ip   = UN[2:2+ni, sj]              # UN[i+1,j]
    un_im   = UN[0:ni, sj]                # UN[i-1,j]
    un_jp   = UN[si, 2:2+nj]              # UN[i,j+1]
    un_jm   = UN[si, 0:nj]                # UN[i,j-1]
    vn_ij   = VN[si, sj]
    vn_ip   = VN[2:2+ni, sj]
    vn_jm   = VN[si, 0:nj]
    vn_ip_jm = VN[2:2+ni, 0:nj]

    f_ij    = F[si, sj]
    f_ip    = F[2:2+ni, sj]

    # --- Advection FUX (x-direction) ---
    sgu = np.where(un_ij >= 0.0, 1.0, -1.0)

    dudr = (un_ip - un_ij) * rx_ip
    dudl = (un_ij - un_im) * rx_i

    rdxa = 1.0 / (dx_i + dx_ip + alpha * sgu * (dx_ip - dx_i))

    fux = rdxa * un_ij * (
        dx_i * dudr + dx_ip * dudl +
        alpha * sgu * (dx_ip * dudl - dx_i * dudr)
    )

    # --- Advection FUY (y-direction) ---
    vbt = (dx_i * vn_ip + dx_ip * vn_ij) * rdelx
    vbb = (dx_i * vn_ip_jm + dx_ip * vn_jm) * rdelx
    vav = 0.5 * (vbt + vbb)

    dyt = 0.5 * (dy_j + dy_jp)
    dyb = 0.5 * (dy_jm + dy_j)

    dudt = (un_jp - un_ij) / dyt
    dudb = (un_ij - un_jm) / dyb

    sgv = np.where(vav >= 0.0, 1.0, -1.0)
    dya = dyt + dyb + alpha * sgv * (dyt - dyb)

    fuy = (vav / dya) * (
        dyb * dudt + dyt * dudb +
        alpha * sgv * (dyt * dudb - dyb * dudt)
    )

    # --- Viscosity (U) ---
    dudxsq = 2.0 * (
        un_im * rx_i / (dx_i + dx_ip) +
        un_ip * rx_ip / (dx_i + dx_ip) -
        un_ij * rx_i * rx_ip
    )

    ubdyt = (dy_j * un_jp + dy_jp * un_ij) / (dy_j + dy_jp)
    ubdyb = (dy_jm * un_ij + dy_j * un_jm) / (dy_j + dy_jm)

    dudyt = (
        un_jp * dy_j * ry_jp -
        un_ij * dy_jp * ry_j -
        ubdyt * (dy_j * ry_jp - dy_jp * ry_j)
    ) / (0.5 * (dy_j + dy_jp))

    dudyb = (
        un_ij * dy_jm * ry_j -
        un_jm * dy_j * ry_jm -
        ubdyb * (dy_jm * ry_j - dy_j * ry_jm)
    ) / (0.5 * (dy_jm + dy_j))

    dudysq = (dudyt - dudyb) * ry_j

    rnucal = f_ij * vnu + (1.0 - f_ij) * vnuc
    visx = rnucal * (dudxsq + dudysq)

    # --- Pressure gradient (U) ---
    p_ij = P[si, sj]
    p_ip = P[2:2+ni, sj]
    rhox = (rhofc + rhod * f_ij) * dx_ip + (rhofc + rhod * f_ip) * dx_i
    gradp_x = np.where(rhox != 0.0, (p_ij - p_ip) * 2.0 / rhox, 0.0)

    # --- Assemble U ---
    fields.U[si, sj] = un_ij + dt * (gradp_x + gx - fux - fuy + visx)

    # ======== V momentum ========
    # Same interior range: i=1..im1-2, j=1..jm1-2
    dx_i_v  = delx[si, np.newaxis]            # delx[i]
    dx_ip_v = delx[2:im1, np.newaxis]         # delx[i+1]
    dx_im_v = delx[0:ni, np.newaxis]          # delx[i-1]
    rx_i_v  = rdx[si, np.newaxis]
    rx_ip_v = rdx[2:im1, np.newaxis]
    rx_im_v = rdx[0:ni, np.newaxis]

    dy_j_v  = dely[np.newaxis, sj]            # dely[j]
    dy_jp_v = dely[np.newaxis, 2:jm1]         # dely[j+1]
    ry_j_v  = rdy[np.newaxis, sj]
    ry_jp_v = rdy[np.newaxis, 2:jm1]

    rdely_v = 1.0 / (dy_j_v + dy_jp_v)

    vn_ij   = VN[si, sj]
    vn_ip   = VN[2:2+ni, sj]
    vn_im   = VN[0:ni, sj]
    vn_jp   = VN[si, 2:2+nj]
    vn_jm   = VN[si, 0:nj]

    un_ij_v  = UN[si, sj]
    un_im_v  = UN[0:ni, sj]
    un_jp_v  = UN[si, 2:2+nj]
    un_im_jp = UN[0:ni, 2:2+nj]

    f_ij_v  = F[si, sj]
    f_jp_v  = F[si, 2:2+nj]

    # --- Advection FVX (x-direction) ---
    ubr = (dy_jp_v * un_ij_v + dy_j_v * un_jp_v) * rdely_v
    ubl = (dy_jp_v * un_im_v + dy_j_v * un_im_jp) * rdely_v
    uav = 0.5 * (ubr + ubl)

    dxr = 0.5 * (dx_i_v + dx_ip_v)
    dxl = 0.5 * (dx_i_v + dx_im_v)

    sgu_v = np.where(uav >= 0.0, 1.0, -1.0)
    dxa = dxr + dxl + alpha * sgu_v * (dxr - dxl)

    dvdr = (vn_ip - vn_ij) / dxr
    dvdl = (vn_ij - vn_im) / dxl

    fvx = (uav / dxa) * (
        dxl * dvdr + dxr * dvdl +
        alpha * sgu_v * (dxr * dvdl - dxl * dvdr)
    )

    # --- Advection FVY (y-direction) ---
    sgv_v = np.where(vn_ij >= 0.0, 1.0, -1.0)
    dya_v = dy_jp_v + dy_j_v + alpha * sgv_v * (dy_jp_v - dy_j_v)

    dvdt = (vn_jp - vn_ij) * ry_jp_v
    dvdb = (vn_ij - vn_jm) * ry_j_v

    fvy = (vn_ij / dya_v) * (
        dy_j_v * dvdt + dy_jp_v * dvdb +
        alpha * sgv_v * (dy_jp_v * dvdb - dy_j_v * dvdt)
    )

    # --- Viscosity (V) ---
    vbdyr = (dx_ip_v * vn_ij + dx_i_v * vn_ip) / (dx_i_v + dx_ip_v)
    vbdyl = (dx_i_v * vn_im + dx_im_v * vn_ij) / (dx_i_v + dx_im_v)

    dvdxr = (
        vn_ip * dx_i_v * rx_ip_v -
        vn_ij * dx_ip_v * rx_i_v -
        vbdyr * (dx_i_v * rx_ip_v - dx_ip_v * rx_i_v)
    ) / (0.5 * (dx_ip_v + dx_i_v))

    dvdxl = (
        vn_ij * dx_im_v * rx_i_v -
        vn_im * dx_i_v * rx_im_v -
        vbdyl * (dx_im_v * rx_i_v - dx_i_v * rx_im_v)
    ) / (0.5 * (dx_i_v + dx_im_v))

    dvdxsq = (dvdxr - dvdxl) * rx_i_v

    dvdysq = 2.0 * (
        vn_jm * ry_j_v / (dy_jp_v + dy_j_v) -
        vn_ij * ry_jp_v * ry_j_v +
        vn_jp * ry_jp_v / (dy_jp_v + dy_j_v)
    )

    rnucal_v = f_ij_v * vnu + (1.0 - f_ij_v) * vnuc
    visy = rnucal_v * (dvdxsq + dvdysq)

    # --- Pressure gradient (V) ---
    p_ij_v = P[si, sj]
    p_jp_v = P[si, 2:2+nj]
    rhoy = (rhofc + rhod * f_ij_v) * dy_jp_v + (rhofc + rhod * f_jp_v) * dy_j_v
    gradp_y = np.where(rhoy != 0.0, (p_ij_v - p_jp_v) * 2.0 / rhoy, 0.0)

    # --- Assemble V ---
    fields.V[si, sj] = vn_ij + dt * (gradp_y + gy - fvx - fvy + visy)

    state.iter = 0


@njit(cache=True)
def _pressit_sweep(U, V, P, PN, F, BETA, delx, dely, rdx, rdy,
                   i_lo, i_hi, j_lo, j_hi,
                   dt, epsi, rhofc, rhod, rcsq, comg):
    """
    Single SOLA pressure sweep (Numba-accelerated).
    Returns True if any cell has |divergence| >= epsi (not yet converged).
    """
    ndelx = delx.shape[0]
    ndely = dely.shape[0]
    not_converged = False

    for j in range(j_lo, j_hi):
        rdy_j = rdy[j]
        dy_j = dely[j]
        dy_jm = dely[j - 1] if j > 0 else dely[0]
        dy_jp = dely[j + 1] if j + 1 < ndely else dely[ndely - 1]

        for i in range(i_lo, i_hi):
            if BETA[i, j] < 0.0:
                continue

            rdx_i = rdx[i]
            dx_i = delx[i]
            dx_im = delx[i - 1] if i > 0 else delx[0]
            dx_ip = delx[i + 1] if i + 1 < ndelx else delx[ndelx - 1]

            # Divergence
            dij = rdx_i * (U[i, j] - U[i - 1, j]) + rdy_j * (V[i, j] - V[i, j - 1])

            # DFUN (compressibility term is zero when rcsq == 0)
            if rcsq != 0.0:
                rhor = rhofc + rhod * F[i, j]
                dfun = dij + rcsq * rhor * (P[i, j] - PN[i, j]) / dt
            else:
                dfun = dij

            # Convergence check
            if abs(dfun) >= epsi:
                not_converged = True

            # Pressure correction (PETA=1, BETA=1 for fluid)
            delp = -BETA[i, j] * dfun
            P[i, j] += delp

            # Velocity corrections: DPTC = 2*DELT*DELP*COMG
            dptc = 2.0 * dt * delp * comg

            # Right: U(i,j)
            i1 = min(i + 1, i_hi)
            if BETA[i1, j] >= 0.0:
                rhoxr = (rhofc + rhod * F[i, j]) * dx_ip + (rhofc + rhod * F[i1, j]) * dx_i
                if rhoxr != 0.0:
                    U[i, j] += dptc / rhoxr

            # Left: U(i-1,j)
            im = max(i - 1, 0)
            if BETA[im, j] >= 0.0:
                rhoxl = (rhofc + rhod * F[im, j]) * dx_i + (rhofc + rhod * F[i, j]) * dx_im
                if rhoxl != 0.0:
                    U[i - 1, j] -= dptc / rhoxl

            # Top: V(i,j)
            j1 = min(j + 1, j_hi)
            if BETA[i, j1] >= 0.0:
                rhoyt = (rhofc + rhod * F[i, j]) * dy_jp + (rhofc + rhod * F[i, j1]) * dy_j
                if rhoyt != 0.0:
                    V[i, j] += dptc / rhoyt

            # Bottom: V(i,j-1)
            jm = max(j - 1, 0)
            if BETA[i, jm] >= 0.0:
                rhoyb = (rhofc + rhod * F[i, jm]) * dy_j + (rhofc + rhod * F[i, j]) * dy_jm
                if rhoyb != 0.0:
                    V[i, j - 1] -= dptc / rhoyb

    return not_converged


def pressure_iteration(mesh: Mesh, fields: Fields, state: RunState,
                       bc_func=None, bc_args=None) -> None:
    """
    SOLA pressure iteration matching Fortran PRESSIT.
    Inner sweep is Numba-JIT compiled for near-Fortran speed.
    """
    U, V, P = fields.U, fields.V, fields.P
    PN = fields.PN
    F = fields.F
    BETA = fields.BETA
    delx, dely = mesh.delx, mesh.dely
    rdx, rdy = mesh.rdx, mesh.rdy

    im1, jm1 = mesh.im1, mesh.jm1
    dt = state.delt
    epsi = state.epsi

    rhofc = state.rhofc
    rhod = state.rhof - state.rhofc

    csq = state.csq
    rcsq = 1.0 / (state.rhof * csq) if csq > 0.0 else 0.0
    comg = 1.0

    state.flg = 1.0
    state.iter = 0
    itmax = 1000

    while state.flg != 0.0:
        if state.iter >= itmax:
            state.fnoc = 1.0
            print(f"PRESSIT: no convergence after {itmax} iterations")
            break
        state.iter += 1

        not_converged = _pressit_sweep(
            U, V, P, PN, F, BETA, delx, dely, rdx, rdy,
            1, im1, 1, jm1,
            dt, epsi, rhofc, rhod, rcsq, comg)

        state.flg = 1.0 if not_converged else 0.0

        if bc_func is not None:
            bc_func(*bc_args)
