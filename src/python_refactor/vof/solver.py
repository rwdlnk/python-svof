# solver.py
import math
import numpy as np
from numba import njit
from .mesh import Mesh
from .fields import Fields
from .state import RunState
from .config import RunConfig

EMF = 1.0e-6
EMF1 = 1.0 - EMF
EM10 = 1.0e-10


@njit(cache=True)
def _tilde_kernel(U, V, UN, VN, P, F, BETA,
                  delx, dely, rdx, rdy, rx, rxi,
                  i_lo, i_hi, j_lo, j_hi,
                  dt, gx, gy, vnu, vnuc,
                  rhof, rhofc, alpha, cyl, nmat):
    """
    Numba kernel for SUBROUTINE TILDE.
    Computes provisional U,V from advection, pressure gradient,
    viscosity, and gravity using the WUDS (alpha) scheme.
    """
    emf = 1.0e-6
    rhod = rhof - rhofc

    for i in range(i_lo, i_hi):
        dx_i  = delx[i]
        dx_ip = delx[i + 1]
        dx_im = delx[i - 1]
        rx_i  = rdx[i]
        rx_ip = rdx[i + 1]
        rx_im = rdx[i - 1]
        rdelx = 1.0 / (dx_i + dx_ip)
        rx_face = rx[i]
        rxi_i = rxi[i]

        for j in range(j_lo, j_hi):
            dy_j  = dely[j]
            dy_jp = dely[j + 1]
            dy_jm = dely[j - 1]
            ry_j  = rdy[j]
            ry_jp = rdy[j + 1]
            ry_jm = rdy[j - 1]

            f_ij = F[i, j]

            # ======== U momentum ========
            f_ip = F[i + 1, j]

            un_ij = UN[i, j]
            un_ip = UN[i + 1, j]
            un_im = UN[i - 1, j]
            un_jp = UN[i, j + 1]
            un_jm = UN[i, j - 1]
            vn_ij = VN[i, j]
            vn_ip = VN[i + 1, j]
            vn_jm_u = VN[i, j - 1]
            vn_ip_jm = VN[i + 1, j - 1]

            # Advection FUX (x-direction)
            sgu = 1.0 if un_ij >= 0.0 else -1.0
            dudr = (un_ip - un_ij) * rx_ip
            dudl = (un_ij - un_im) * rx_i
            rdxa = 1.0 / (dx_i + dx_ip + alpha * sgu * (dx_ip - dx_i))
            fux = rdxa * un_ij * (
                dx_i * dudr + dx_ip * dudl +
                alpha * sgu * (dx_ip * dudl - dx_i * dudr))

            # Advection FUY (y-direction)
            vbt = (dx_i * vn_ip + dx_ip * vn_ij) * rdelx
            vbb = (dx_i * vn_ip_jm + dx_ip * vn_jm_u) * rdelx
            vav = 0.5 * (vbt + vbb)
            dyt = 0.5 * (dy_j + dy_jp)
            dyb = 0.5 * (dy_jm + dy_j)
            dudt_u = (un_jp - un_ij) / dyt
            dudb_u = (un_ij - un_jm) / dyb
            sgv = 1.0 if vav >= 0.0 else -1.0
            dya = dyt + dyb + alpha * sgv * (dyt - dyb)
            fuy = (vav / dya) * (
                dyb * dudt_u + dyt * dudb_u +
                alpha * sgv * (dyt * dudb_u - dyb * dudt_u))

            # Viscosity (U)
            dudxsq = 2.0 * (
                un_im * rx_i * rdelx +
                un_ip * rx_ip * rdelx -
                un_ij * rx_i * rx_ip)

            dy_j_jp = dy_j + dy_jp
            dy_jm_j = dy_jm + dy_j
            ubdyt = (dy_j * un_jp + dy_jp * un_ij) / dy_j_jp
            ubdyb = (dy_jm * un_ij + dy_j * un_jm) / dy_jm_j

            dudyt = (un_jp * dy_j * ry_jp -
                     un_ij * dy_jp * ry_j -
                     ubdyt * (dy_j * ry_jp - dy_jp * ry_j)) / (0.5 * dy_j_jp)

            dudyb = (un_ij * dy_jm * ry_j -
                     un_jm * dy_j * ry_jm -
                     ubdyb * (dy_jm * ry_j - dy_j * ry_jm)) / (0.5 * dy_jm_j)

            dudysq = (dudyt - dudyb) * ry_j

            rnucal = f_ij * vnu + (1.0 - f_ij) * vnuc
            visx = rnucal * (dudxsq + dudysq)

            # CYL viscous terms for U
            if cyl > 0.5:
                dudxl_c = (un_ij - un_im) * rx_i
                dudxr_c = (un_ip - un_ij) * rx_ip
                rxdudx = rx_face * (dx_ip * dudxl_c + dx_i * dudxr_c) * rdelx
                rxsqu = un_ij * rx_face * rx_face
                visx = rnucal * (dudxsq + dudysq + cyl * rxdudx - cyl * rxsqu)

            # Pressure gradient (U)
            p_ij = P[i, j]
            p_ip = P[i + 1, j]
            rhox = (rhofc + rhod * f_ij) * dx_ip + (rhofc + rhod * f_ip) * dx_i
            gradp_x = (p_ij - p_ip) * 2.0 / rhox if rhox != 0.0 else 0.0

            # Assemble U
            u_new = un_ij + dt * (gradp_x + gx - fux - fuy + visx)
            beta_ok = (BETA[i, j] >= 0.0) and (BETA[i + 1, j] >= 0.0)
            if nmat == 1:
                beta_ok = beta_ok and ((f_ij + f_ip) >= emf)
            U[i, j] = u_new if beta_ok else 0.0

            # ======== V momentum ========
            f_jp = F[i, j + 1]

            vn_ij_v = VN[i, j]
            vn_ip_v = VN[i + 1, j]
            vn_im_v = VN[i - 1, j]
            vn_jp_v = VN[i, j + 1]
            vn_jm_v = VN[i, j - 1]
            un_ij_v = UN[i, j]
            un_im_v = UN[i - 1, j]
            un_jp_v = UN[i, j + 1]
            un_im_jp = UN[i - 1, j + 1]

            rdely = 1.0 / (dy_j + dy_jp)

            # Advection FVX (x-direction)
            ubr = (dy_jp * un_ij_v + dy_j * un_jp_v) * rdely
            ubl = (dy_jp * un_im_v + dy_j * un_im_jp) * rdely
            uav = 0.5 * (ubr + ubl)
            dxr = 0.5 * (dx_i + dx_ip)
            dxl = 0.5 * (dx_i + dx_im)
            sgu_v = 1.0 if uav >= 0.0 else -1.0
            dxa_v = dxr + dxl + alpha * sgu_v * (dxr - dxl)
            dvdr = (vn_ip_v - vn_ij_v) / dxr
            dvdl = (vn_ij_v - vn_im_v) / dxl
            fvx = (uav / dxa_v) * (
                dxl * dvdr + dxr * dvdl +
                alpha * sgu_v * (dxr * dvdl - dxl * dvdr))

            # Advection FVY (y-direction)
            sgv_v = 1.0 if vn_ij_v >= 0.0 else -1.0
            dya_v = dy_jp + dy_j + alpha * sgv_v * (dy_jp - dy_j)
            dvdt_v = (vn_jp_v - vn_ij_v) * ry_jp
            dvdb_v = (vn_ij_v - vn_jm_v) * ry_j
            fvy = (vn_ij_v / dya_v) * (
                dy_j * dvdt_v + dy_jp * dvdb_v +
                alpha * sgv_v * (dy_jp * dvdb_v - dy_j * dvdt_v))

            # Viscosity (V)
            dx_i_ip = dx_i + dx_ip
            dx_i_im = dx_i + dx_im
            vbdyr = (dx_ip * vn_ij_v + dx_i * vn_ip_v) / dx_i_ip
            vbdyl = (dx_i * vn_im_v + dx_im * vn_ij_v) / dx_i_im

            dvdxr_v = (vn_ip_v * dx_i * rx_ip -
                       vn_ij_v * dx_ip * rx_i -
                       vbdyr * (dx_i * rx_ip - dx_ip * rx_i)) / (0.5 * dx_i_ip)

            dvdxl_v = (vn_ij_v * dx_im * rx_i -
                       vn_im_v * dx_i * rx_im -
                       vbdyl * (dx_im * rx_i - dx_i * rx_im)) / (0.5 * dx_i_im)

            dvdxsq = (dvdxr_v - dvdxl_v) * rx_i

            dy_jp_j = dy_jp + dy_j
            dvdysq = 2.0 * (
                vn_jm_v * ry_j / dy_jp_j -
                vn_ij_v * ry_jp * ry_j +
                vn_jp_v * ry_jp / dy_jp_j)

            rnucal_v = f_ij * vnu + (1.0 - f_ij) * vnuc
            visy = rnucal_v * (dvdxsq + dvdysq)

            # CYL viscous terms for V
            if cyl > 0.5:
                dvdxrx = (vbdyr - vbdyl) * rx_i * rxi_i
                visy = rnucal_v * (dvdxsq + dvdysq + cyl * dvdxrx)

            # Pressure gradient (V)
            p_jp = P[i, j + 1]
            rhoy = (rhofc + rhod * f_ij) * dy_jp + (rhofc + rhod * f_jp) * dy_j
            gradp_y = (p_ij - p_jp) * 2.0 / rhoy if rhoy != 0.0 else 0.0

            # Assemble V
            v_new = vn_ij_v + dt * (gradp_y + gy - fvx - fvy + visy)
            beta_ok_v = (BETA[i, j] >= 0.0) and (BETA[i, j + 1] >= 0.0)
            if nmat == 1:
                beta_ok_v = beta_ok_v and ((f_ij + f_jp) >= emf)
            V[i, j] = v_new if beta_ok_v else 0.0


def tilde_step(mesh: Mesh, fields: Fields, state: RunState) -> None:
    """
    Numba-accelerated version of Fortran SUBROUTINE TILDE for U and V.
    Computes provisional velocities from advection, pressure gradient,
    viscosity, and gravity using the WUDS (alpha) scheme.
    """
    _tilde_kernel(
        fields.U, fields.V, fields.UN, fields.VN,
        fields.P, fields.F, fields.BETA,
        mesh.delx, mesh.dely, mesh.rdx, mesh.rdy, mesh.rx, mesh.rxi,
        1, mesh.im1, 1, mesh.jm1,
        state.delt, state.gx, state.gy, state.vnu, state.vnuc,
        state.rhof, state.rhofc, state.alpha, state.cyl, state.nmat)

    state.iter = 0


@njit(cache=True)
def _pressit_sweep(U, V, P, PN, F, BETA, PETA, NF, PS, PR,
                   delx, dely, rdx, rdy, rxi,
                   i_lo, i_hi, j_lo, j_hi,
                   dt, epsi, rhofc, rhod, rcsq, comg, nmat, cyl):
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

            # NMAT=1: skip empty cells
            if nmat != 2 and F[i, j] < EMF:
                continue

            rdx_i = rdx[i]
            dx_i = delx[i]
            dx_im = delx[i - 1] if i > 0 else delx[0]
            dx_ip = delx[i + 1] if i + 1 < ndelx else delx[ndelx - 1]

            # --- NMAT=1 surface cell pressure path ---
            nf_ij = NF[i, j]
            if nmat != 2 and nf_ij != 0 and nf_ij != 5:
                # Surface cell: compute pressure from interpolation
                l_idx = i
                m_idx = j
                if nf_ij == 1:
                    l_idx = i - 1
                elif nf_ij == 2:
                    l_idx = i + 1
                elif nf_ij == 3:
                    m_idx = j - 1
                elif nf_ij == 4:
                    m_idx = j + 1

                # Compute NFE = max of neighbor NF values
                nfel = NF[i - 1, j]
                nfer = NF[i + 1, j] if i + 1 < NF.shape[0] else 0
                nfeb = NF[i, j - 1]
                nfet = NF[i, j + 1] if j + 1 < NF.shape[1] else 0
                nfe = max(nfel, nfer, nfeb, nfet)
                nfe = min(nfe, len(PR) - 1)

                psurf = PS[i, j] + PR[nfe]
                plm = P[l_idx, m_idx]
                if NF[l_idx, m_idx] != 0 and BETA[i, j] > 0.0:
                    plm = psurf
                delp = (1.0 - PETA[i, j]) * plm + PETA[i, j] * psurf - P[i, j]
            else:
                # --- Normal divergence-based path ---
                # Divergence with CYL term
                nrxi = rxi.shape[0]
                dij = rdx_i * (U[i, j] - U[i - 1, j]) + rdy_j * (V[i, j] - V[i, j - 1])
                if cyl > 0.5 and i < nrxi:
                    dij += cyl * rxi[i] * (U[i, j] + U[i - 1, j]) / 2.0

                # DFUN (compressibility term is zero when rcsq == 0)
                if rcsq != 0.0:
                    rhof_local = rhofc + rhod * F[i, j]
                    rhor = (rhofc + rhod) / rhof_local  # RHOF/rho(i,j)
                    dfun = dij + rhor * rcsq * (P[i, j] - PN[i, j]) / dt
                else:
                    dfun = dij

                # Convergence check
                if abs(dfun) >= epsi:
                    not_converged = True

                # Pressure correction: DELP = -BETA * DFUN * PETA
                delp = -BETA[i, j] * dfun * PETA[i, j]

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

    # Compute COMG from RDTEXP (matches Fortran PRESSIT lines 1904-1905)
    ctos = dt * state.rdtexp
    comg = min(ctos * ctos, 1.0)

    PETA = fields.PETA

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
            U, V, P, PN, F, BETA, PETA, fields.NF, fields.PS, fields.PR,
            delx, dely, rdx, rdy, mesh.rxi,
            1, im1, 1, jm1,
            dt, epsi, rhofc, rhod, rcsq, comg, state.nmat, state.cyl)

        state.flg = 1.0 if not_converged else 0.0

        if bc_func is not None:
            bc_func(*bc_args)


# =====================================================================
# VOF advection (VFCONV) — Hirt-Nichols donor-acceptor method
# =====================================================================

@njit(cache=True)
def _vfconv_kernel(U, V, F, FN, NF, BETA, delx, dely, rdx, rdy, xi,
                   im1, jm1, imax, jmax, dt, cyl, i_lo, j_lo,
                   skip_phase_b):
    """
    Numba kernel for VFCONV.  Operates on 0-based arrays.
    Fortran indices 1..IM1, 1..JM1  →  Python 0..im1-1, 0..jm1-1
    Fortran IMAX (1-based size)     →  Python imax (0-based size)

    i_lo, j_lo: starting indices for Phase A loops.  Use 1 to skip ghost
    cells (required for MPI to avoid double-counting flux at inter-rank
    boundaries; safe for serial since ghost-cell fluxes are zero at walls
    and their F values are overwritten by BC anyway).

    skip_phase_b: if True, skip Phase B (F cleanup/clipping).  Used by
    MPI to defer Phase B until after halo_accumulate + halo_exchange.

    Returns (flgc, vchgt_delta).
    """
    PI = 3.14159265358979
    flgc = 0.0
    vchgt = 0.0

    # --- Phase A: donor-acceptor flux ---
    for j in range(j_lo, jm1):
        for i in range(i_lo, im1):
            vx = U[i, j] * dt
            vy = V[i, j] * dt
            abvx = abs(vx)
            abvy = abs(vy)

            if abvx > 0.5 * delx[i] or abvy > 0.5 * dely[j]:
                flgc = 1.0

            # --- X-direction flux ---
            nxi = xi.shape[0]
            if vx >= 0.0:
                ia = i + 1
                id_ = i
                idm = max(i - 1, 0)
                rb = xi[i] if i < nxi else xi[nxi - 1]  # x[i] in Fortran = boundary face
                ra = xi[min(i + 1, nxi - 1)]
                rd = xi[i]
            else:
                ia = i
                id_ = i + 1
                idm = min(i + 2, imax - 1)
                rb = xi[i] if i < nxi else xi[nxi - 1]
                ra = xi[i]
                rd = xi[min(i + 1, nxi - 1)]

            iad = ia
            if NF[id_, j] == 3 or NF[id_, j] == 4:
                iad = id_
            if FN[ia, j] < EMF or FN[idm, j] < EMF:
                iad = ia

            fdm = max(FN[idm, j], FN[id_, j])
            fx1 = FN[iad, j] * abvx + max(
                (fdm - FN[iad, j]) * abvx - (fdm - FN[id_, j]) * delx[id_],
                0.0)
            fx = min(fx1, FN[id_, j] * delx[id_])

            # CYL terms: radius-ratio flux weighting
            if cyl > 0.5 and abs(rd) > 1.0e-30:
                ratio_d = abs(rb / rd) * cyl + (1.0 - cyl)
                ratio_a = abs(rb / ra) * cyl + (1.0 - cyl) if abs(ra) > 1.0e-30 else 1.0
            else:
                ratio_d = 1.0
                ratio_a = 1.0

            F[id_, j] -= fx * rdx[id_] * ratio_d
            F[ia, j] += fx * rdx[ia] * ratio_a

            # --- Y-direction flux ---
            if vy >= 0.0:
                ja = j + 1
                jd = j
                jdm = max(j - 1, 0)
            else:
                ja = j
                jd = j + 1
                jdm = min(j + 2, jmax - 1)

            jad = ja
            if NF[i, jd] == 1 or NF[i, jd] == 2:
                jad = jd
            if FN[i, ja] < EMF or FN[i, jdm] < EMF:
                jad = ja

            fdm = max(FN[i, jdm], FN[i, jd])
            fy1 = FN[i, jad] * abvy + max(
                (fdm - FN[i, jad]) * abvy - (fdm - FN[i, jd]) * dely[jd],
                0.0)
            fy = min(fy1, FN[i, jd] * dely[jd])

            F[i, jd] -= fy * rdy[jd]
            F[i, ja] += fy * rdy[ja]

    # --- Phase B: F cleanup ---
    if skip_phase_b:
        return flgc, vchgt

    # Fortran loops I=2..IM1, J=2..JM1 (1-based) → 1..im1-1, 1..jm1-1 (0-based)
    for j in range(1, jm1):
        for i in range(1, im1):
            if BETA[i, j] < 0.0:
                continue

            vchg = 0.0
            if F[i, j] > EMF and F[i, j] < EMF1:
                pass  # surface cell, no snapping
            elif F[i, j] >= EMF1:
                vchg = -(1.0 - F[i, j])
                F[i, j] = 1.0
            else:
                vchg = F[i, j]
                F[i, j] = 0.0

            # Volume element: planar or cylindrical
            vol_factor = xi[i] * 2.0 * PI * cyl + (1.0 - cyl) if i < nxi else (1.0 - cyl)
            vchgt += vchg * delx[i] * dely[j] * vol_factor

            # Full cell adjacent to void: bleed off to maintain surface layer
            if F[i, j] >= EMF1:
                if (F[i + 1, j] < EMF or F[i - 1, j] < EMF or
                        F[i, j + 1] < EMF or F[i, j - 1] < EMF):
                    F[i, j] -= 1.1 * EMF
                    vchg2 = 1.1 * EMF
                    vchgt += vchg2 * delx[i] * dely[j] * vol_factor

    return flgc, vchgt


def vfconv(mesh: Mesh, fields: Fields, state: RunState,
           i_lo: int = 1, j_lo: int = 1,
           skip_phase_b: bool = False) -> None:
    """
    VOF advection: port of Fortran SUBROUTINE VFCONV.
    Convects the F field using the Hirt-Nichols donor-acceptor method.

    i_lo, j_lo: starting indices for Phase A loops (default 1, skipping
    ghost cells — safe for serial and required for MPI).

    skip_phase_b: if True, only run Phase A (flux computation).  Phase B
    (F cleanup/clipping) must be called separately via vfconv_phase_b().
    Used by MPI to insert halo_accumulate+halo_exchange between phases.
    """
    if state.icyle < 1:
        return

    flgc, vchgt_delta = _vfconv_kernel(
        fields.U, fields.V, fields.F, fields.FN, fields.NF, fields.BETA,
        mesh.delx, mesh.dely, mesh.rdx, mesh.rdy, mesh.xi,
        mesh.im1, mesh.jm1, mesh.imax, mesh.jmax,
        state.delt, state.cyl, i_lo, j_lo, skip_phase_b)

    state.flgc = flgc
    state.vchgt += vchgt_delta
    if flgc > 0.5:
        state.nflgc += 1


@njit(cache=True)
def _vfconv_phase_b_kernel(F, BETA, delx, dely, xi, im1, jm1, cyl):
    """Phase B of vfconv: F cleanup (clipping and bleedoff)."""
    PI = 3.14159265358979
    nxi = xi.shape[0]
    vchgt = 0.0

    for j in range(1, jm1):
        for i in range(1, im1):
            if BETA[i, j] < 0.0:
                continue

            vchg = 0.0
            if F[i, j] > EMF and F[i, j] < EMF1:
                pass  # surface cell, no snapping
            elif F[i, j] >= EMF1:
                vchg = -(1.0 - F[i, j])
                F[i, j] = 1.0
            else:
                vchg = F[i, j]
                F[i, j] = 0.0

            vol_factor = xi[i] * 2.0 * PI * cyl + (1.0 - cyl) if i < nxi else (1.0 - cyl)
            vchgt += vchg * delx[i] * dely[j] * vol_factor

            # Full cell adjacent to void: bleed off to maintain surface layer
            if F[i, j] >= EMF1:
                if (F[i + 1, j] < EMF or F[i - 1, j] < EMF or
                        F[i, j + 1] < EMF or F[i, j - 1] < EMF):
                    F[i, j] -= 1.1 * EMF
                    vchg2 = 1.1 * EMF
                    vchgt += vchg2 * delx[i] * dely[j] * vol_factor

    return vchgt


def vfconv_phase_b(mesh: Mesh, fields: Fields, state: RunState) -> None:
    """Run Phase B of vfconv (F cleanup) as a separate step.

    Called after halo_accumulate + halo_exchange in MPI mode so that
    F clipping sees the complete flux from all neighbors.
    """
    vchgt_delta = _vfconv_phase_b_kernel(
        fields.F, fields.BETA,
        mesh.delx, mesh.dely, mesh.xi,
        mesh.im1, mesh.jm1, state.cyl)
    state.vchgt += vchgt_delta


# =====================================================================
# PETACAL — compute NF flags (simplified for NMAT=2)
# =====================================================================

@njit(cache=True)
def _petacal_kernel(F, P, NF, PS, PETA, BETA, DTANTH,
                    delx, dely, rdx, rdy, xi, x, rxi,
                    im1, jm1, imax, jmax,
                    iwl, iwr, iwt, iwb,
                    isurf10, sigma, nmat, cyl):
    """
    Numba kernel for PETACAL.  Computes NF flags, curvature, surface
    pressure PS, and PETA interpolation factor.
    Fortran loops I=2..IM1, J=2..JM1 (1-based) → 1..im1-1, 1..jm1-1 (0-based).
    """
    EP10 = 1.0e10

    # Initialize all arrays
    for j in range(jmax):
        for i in range(imax):
            NF[i, j] = 0
            PS[i, j] = 0.0
            PETA[i, j] = 1.0

    for j in range(1, jm1):
        for i in range(1, im1):
            DTANTH[i, j] = EP10

            if BETA[i, j] < 0.0:
                continue

            if F[i, j] < EMF:
                NF[i, j] = 6
            if F[i, j] < EMF or F[i, j] > EMF1:
                continue

            # Check if this is a surface cell (has an empty neighbor)
            if (F[i + 1, j] >= EMF and F[i, j + 1] >= EMF and
                    F[i - 1, j] >= EMF and F[i, j - 1] >= EMF):
                continue

            # --- This is a surface cell: compute F gradients ---
            dxr = (delx[i] + delx[i + 1]) / 2.0
            dxl = (delx[i] + delx[i - 1]) / 2.0
            dyt = (dely[j] + dely[j + 1]) / 2.0
            dyb = (dely[j] + dely[j - 1]) / 2.0

            rxden = 1.0 / (dxr * dxl * (dxr + dxl))
            ryden = 1.0 / (dyt * dyb * (dyt + dyb))

            # Area-weighted F averages from 3x3 stencil
            # Top row (j+1)
            fl = F[i - 1, j + 1]
            if BETA[i - 1, j + 1] < 0.0 or (i == 1 and iwl < 3):
                fl = 1.0
            fc = F[i, j + 1]
            if BETA[i, j + 1] < 0.0:
                fc = 1.0
            fr = F[i + 1, j + 1]
            if BETA[i + 1, j + 1] < 0.0 or (i == im1 - 1 and iwr < 3):
                fr = 1.0
            avft = fl * delx[i - 1] + fc * delx[i] + fr * delx[i + 1]

            # Bottom row (j-1)
            fl = F[i - 1, j - 1]
            if BETA[i - 1, j - 1] < 0.0 or (i == 1 and iwl < 3):
                fl = 1.0
            fc = F[i, j - 1]
            if BETA[i, j - 1] < 0.0:
                fc = 1.0
            fr = F[i + 1, j - 1]
            if BETA[i + 1, j - 1] < 0.0 or (i == im1 - 1 and iwr < 3):
                fr = 1.0
            avfb = fl * delx[i - 1] + fc * delx[i] + fr * delx[i + 1]

            # Center row (j) — x-average
            fl = F[i - 1, j]
            if BETA[i - 1, j] < 0.0 or (i == 1 and iwl < 3):
                fl = 1.0
            fr = F[i + 1, j]
            if BETA[i + 1, j] < 0.0 or (i == im1 - 1 and iwr < 3):
                fr = 1.0
            avfcy = fl * delx[i - 1] + F[i, j] * delx[i] + fr * delx[i + 1]

            # Center column (i) — y-average
            fb = F[i, j - 1]
            if BETA[i, j - 1] < 0.0 or (j == 1 and iwb < 3):
                fb = 1.0
            ft = F[i, j + 1]
            if BETA[i, j + 1] < 0.0 or (j == jm1 - 1 and iwt < 3):
                ft = 1.0
            avfcx = fb * dely[j - 1] + F[i, j] * dely[j] + ft * dely[j + 1]

            # Left column (i-1) — y-average
            fb = F[i - 1, j - 1]
            if BETA[i - 1, j - 1] < 0.0 or (j == 1 and iwb < 3):
                fb = 1.0
            fc = F[i - 1, j]
            if BETA[i - 1, j] < 0.0:
                fc = 1.0
            ft = F[i - 1, j + 1]
            if BETA[i - 1, j + 1] < 0.0 or (j == jm1 - 1 and iwt < 3):
                ft = 1.0
            avfl = fb * dely[j - 1] + fc * dely[j] + ft * dely[j + 1]

            # Right column (i+1) — y-average
            fb = F[i + 1, j - 1]
            if BETA[i + 1, j - 1] < 0.0 or (j == 1 and iwb < 3):
                fb = 1.0
            fc = F[i + 1, j]
            if BETA[i + 1, j] < 0.0:
                fc = 1.0
            ft = F[i + 1, j + 1]
            if BETA[i + 1, j + 1] < 0.0 or (j == jm1 - 1 and iwt < 3):
                ft = 1.0
            avfr = fb * dely[j - 1] + fc * dely[j] + ft * dely[j + 1]

            # Thickness measures for surface tension validity check
            xthm = 3.0 * max(avft, avfcy, avfb) / (delx[i - 1] + delx[i] + delx[i + 1])
            ythm = 3.0 * max(avfl, avfcx, avfr) / (dely[j - 1] + dely[j] + dely[j + 1])

            # F gradients
            pfx = rxden * ((avfr - avfcx) * dxl * dxl +
                           (avfcx - avfl) * dxr * dxr)
            pfy = ryden * ((avft - avfcy) * dyb * dyb +
                           (avfcy - avfb) * dyt * dyt)

            pf = pfx * pfx + pfy * pfy
            if pf <= EM10:
                NF[i, j] = 5
                P[i, j] = (P[i + 1, j] + P[i, j + 1] +
                           P[i - 1, j] + P[i, j - 1]) / 4.0
                continue

            # Determine NF direction, L, M, DMX, DMIN from gradient
            abpfx = abs(pfx)
            abpfy = abs(pfy)

            # L, M: coordinates of adjacent full cell for pressure interpolation
            l_idx = i
            m_idx = j
            dmx = 0.0
            dmin = 0.0
            dxdyr = 0.0

            if abpfy < abpfx:
                # x-dominant gradient
                dxdyr = dely[j] * rdx[i]
                if pfx > 0.0:
                    NF[i, j] = 2
                    DTANTH[i, j] = -pfy
                    l_idx = i + 1
                    dmx = delx[i]
                    dmin = (dmx + delx[i + 1]) / 2.0
                else:
                    NF[i, j] = 1
                    DTANTH[i, j] = pfy
                    l_idx = i - 1
                    dmx = delx[i]
                    dmin = (dmx + delx[i - 1]) / 2.0
            else:
                # y-dominant gradient
                dxdyr = delx[i] * rdy[j]
                if pfy > 0.0:
                    NF[i, j] = 4
                    DTANTH[i, j] = -pfx
                    m_idx = j + 1
                    dmx = dely[j]
                    dmin = (dmx + dely[j + 1]) / 2.0
                else:
                    NF[i, j] = 3
                    DTANTH[i, j] = pfx
                    m_idx = j - 1
                    dmx = dely[j]
                    dmin = (dmx + dely[j - 1]) / 2.0

            abdtan = abs(DTANTH[i, j])

            # --- Curvature and surface pressure (Fortran lines 1587-1638) ---
            dfs = (0.5 - F[i, j]) * dmx
            if F[i, j] < 0.5 * abdtan * dxdyr:
                arg = 8.0 * F[i, j] * dxdyr * abdtan
                dfs = dmx * (1.0 + dxdyr * abdtan - arg ** 0.5) / 2.0

            if isurf10 >= 1:
                nfc = NF[i, j]
                pxr = (avfr - avfcx) / dxr
                pxl = (avfcx - avfl) / dxl
                pyt = (avft - avfcy) / dyt
                pyb = (avfcy - avfb) / dyb
                ydfs = -dfs
                if nfc == 2 or nfc == 4:
                    ydfs = dfs

                if nfc <= 2:
                    # x-dominant: compute slopes in y-direction
                    dxdn = dely[j]
                    xinb = ydfs + DTANTH[i, j] * dxdn / 2.0
                    xint = 2.0 * ydfs - xinb
                    gp1 = pyt
                    px1 = pxl
                    if xint > 0.0:
                        px1 = pxr
                    if abs(px1) < abs(gp1):
                        gp1 = (1.0 if gp1 >= 0.0 else -1.0) / (abs(px1) + EM10)
                    gp2 = pyb
                    px2 = pxr
                    if xinb < 0.0:
                        px2 = pxl
                    if abs(px2) < abs(gp2):
                        gp2 = (1.0 if gp2 >= 0.0 else -1.0) / (abs(px2) + EM10)
                else:
                    # y-dominant: compute slopes in x-direction
                    dxdn = delx[i]
                    yinr = ydfs + DTANTH[i, j] * dxdn / 2.0
                    yinl = 2.0 * ydfs - yinr
                    gp1 = pxr
                    py1 = pyt
                    if yinr < 0.0:
                        py1 = pyb
                    if abs(py1) < abs(gp1):
                        gp1 = (1.0 if gp1 >= 0.0 else -1.0) / (abs(py1) + EM10)
                    gp2 = pxl
                    py2 = pyb
                    if yinl > 0.0:
                        py2 = pyt
                    if abs(py2) < abs(gp2):
                        gp2 = (1.0 if gp2 >= 0.0 else -1.0) / (abs(py2) + EM10)

                gp1d = 1.0 + gp1 * gp1
                gp2d = 1.0 + gp2 * gp2
                curvxy = (gp2 / gp2d ** 0.5 - gp1 / gp1d ** 0.5) / dxdn
                # CYL curvature
                curvcyl = 0.0
                if cyl >= 1.0:
                    nfc = NF[i, j]
                    nxi = xi.shape[0]
                    nx = x.shape[0]
                    xlitlr = xi[i] if i < nxi else xi[nxi - 1]
                    if nfc == 1 and i - 1 >= 0 and i - 1 < nx:
                        xlitlr = x[i - 1] + F[i, j] * delx[i]
                    elif nfc == 2 and i < nx:
                        xlitlr = x[i] - F[i, j] * delx[i]
                    nrxi = rxi.shape[0]
                    if abs(xlitlr) > EM10:
                        rlitlr = 1.0 / xlitlr
                        if nrxi > 1:
                            rlitlr = min(rlitlr, rxi[1])
                    else:
                        rlitlr = rxi[1] if nrxi > 1 else 0.0
                    abdtan_val = abs(DTANTH[i, j])
                    atan_val = math.atan(abdtan_val)
                    if nfc <= 2:
                        trig = abs(math.cos(atan_val))
                    else:
                        trig = abs(math.sin(atan_val))
                    pfx_sign = 1.0 if pfx >= 0.0 else -1.0
                    curvcyl = -cyl * trig * pfx_sign * rlitlr
                curv = curvxy + curvcyl
                PS[i, j] = sigma * curv
                if xthm < 1.0 or ythm < 1.0:
                    PS[i, j] = 0.0

            # --- NFSB: check if all neighbors are void/obstacle ---
            nfsb = 0
            if F[i + 1, j] < EMF or i == im1 - 1 or BETA[i + 1, j] < 0.0:
                nfsb += 1
            if F[i, j + 1] < EMF or BETA[i, j + 1] < 0.0:
                nfsb += 2
            if F[i - 1, j] < EMF or BETA[i - 1, j] < 0.0:
                nfsb += 4
            if F[i, j - 1] < EMF or BETA[i, j - 1] < 0.0:
                nfsb += 8
            if nfsb == 15:
                PS[i, j] = 0.0

            # --- PETA for NMAT=1 (free surface) ---
            if nmat != 2 and dmin > 0.0:
                dfs_over_dmin = dfs / dmin
                if dfs_over_dmin < 1.0:
                    PETA[i, j] = 1.0 / (1.0 - dfs_over_dmin)
                if l_idx == 0 or l_idx == imax - 1:
                    PETA[i, j] = 1.0
                if m_idx == 0 or m_idx == jmax - 1:
                    PETA[i, j] = 1.0
                if BETA[l_idx, m_idx] < 0.0:
                    PETA[i, j] = 1.0


def petacal(mesh: Mesh, fields: Fields, state: RunState,
            cfg: RunConfig) -> None:
    """
    Compute NF flags, curvature, surface pressure PS, and PETA for
    surface cells.  Port of Fortran PETACAL.
    """
    _petacal_kernel(
        fields.F, fields.P, fields.NF, fields.PS, fields.PETA,
        fields.BETA, fields.DTANTH,
        mesh.delx, mesh.dely, mesh.rdx, mesh.rdy,
        mesh.xi, mesh.x, mesh.rxi,
        mesh.im1, mesh.jm1, mesh.imax, mesh.jmax,
        cfg.iwl, cfg.iwr, cfg.iwt, cfg.iwb,
        1 if cfg.isurf10 else 0, state.sigma, cfg.nmat, state.cyl)


# =====================================================================
# TMS10 — two-material surface tension (port of Fortran TMS10)
# =====================================================================

@njit(cache=True)
def _tms10_kernel(U, V, F, NF, PS, BETA, DTANTH, delx, dely,
                  im1, jm1, dt, rhofc, rhod):
    """
    Numba kernel for TMS10: apply surface tension forces to velocities
    at surface cells for NMAT=2.
    """
    for j in range(1, jm1):
        for i in range(1, im1):
            nf_ij = NF[i, j]
            if nf_ij == 0 or nf_ij >= 5 or BETA[i, j] < 0.0:
                continue

            whtl = 0.5
            whtr = 0.5
            whtt = 0.5
            whtb = 0.5

            if nf_ij <= 2:
                # x-dominant interface
                whtl = 1.0 - F[i, j]
                if nf_ij == 2:
                    whtl = 1.0 - whtl
                whtr = 1.0 - whtl
                stfx = PS[i, j] * dely[j]
                if nf_ij == 1:
                    stfx = -stfx
                stfy = stfx * DTANTH[i, j]
            else:
                # y-dominant interface
                whtb = 1.0 - F[i, j]
                if nf_ij == 4:
                    whtb = 1.0 - whtb
                whtt = 1.0 - whtb
                stfy = PS[i, j] * delx[i]
                if nf_ij == 3:
                    stfy = -stfy
                stfx = -stfy * DTANTH[i, j]

            ndelx = delx.shape[0]
            ndely = dely.shape[0]
            dx_ip = delx[min(i + 1, ndelx - 1)]
            dx_im = delx[max(i - 1, 0)]
            dy_jp = dely[min(j + 1, ndely - 1)]
            dy_jm = dely[max(j - 1, 0)]

            rhoxr = (rhofc + rhod * F[i, j]) * dx_ip + (rhofc + rhod * F[i + 1, j]) * delx[i]
            U[i, j] += 2.0 * dt * whtr * stfx / (rhoxr * dely[j])

            rhoxl = (rhofc + rhod * F[i - 1, j]) * delx[i] + (rhofc + rhod * F[i, j]) * dx_im
            U[i - 1, j] += 2.0 * dt * whtl * stfx / (rhoxl * dely[j])

            rhoyt = (rhofc + rhod * F[i, j]) * dy_jp + (rhofc + rhod * F[i, j + 1]) * dely[j]
            V[i, j] += 2.0 * dt * whtt * stfy / (rhoyt * delx[i])

            rhoyb = (rhofc + rhod * F[i, j - 1]) * dely[j] + (rhofc + rhod * F[i, j]) * dy_jm
            V[i, j - 1] += 2.0 * dt * whtb * stfy / (rhoyb * delx[i])


def tms10(mesh, fields, state):
    """
    Two-material surface tension.  Port of Fortran TMS10.
    Adds surface tension forces to U,V at NMAT=2 interface cells.
    Called between tilde_step and BC (Fortran line 138).
    """
    _tms10_kernel(
        fields.U, fields.V, fields.F, fields.NF, fields.PS,
        fields.BETA, fields.DTANTH,
        mesh.delx, mesh.dely,
        mesh.im1, mesh.jm1, state.delt,
        state.rhofc, state.rhof - state.rhofc)


# =====================================================================
# BETA computation — port of Fortran SETUP lines 2641-2679
# =====================================================================

@njit(cache=True)
def _compute_beta_kernel(F, BETA, delx, dely, rdx, rdy,
                         im1, jm1, delt,
                         rhof, rhofc, rhod, csq, omg, vnu, vnuc, sigma, cyl):
    """
    Numba kernel for computing BETA(i,j) and stability limits.
    Fortran indices 2..IM1 → Python 1..im1-1 (0-based).
    Returns (ds, dtvis, dtsft, rdtexp).
    """
    ds = 1.0e10
    dtvis = 1.0e10

    ndelx = delx.shape[0]
    ndely = dely.shape[0]

    # --- DTVIS and DS ---
    for j in range(1, jm1):
        for i in range(1, im1):
            dxsq = delx[i] * delx[i]
            dysq = dely[j] * dely[j]
            rdsq = dxsq * dysq / (dxsq + dysq)
            rnumax = max(vnu, vnuc)
            rdsq = rdsq / (3.0 * rnumax + 1.0e-10)
            dtvis = min(dtvis, rdsq)
            ds = min(delx[i], dely[j], ds)

    # --- DTSFT ---
    sigx = sigma
    rhomn = min(rhof, rhofc)
    if sigx == 0.0:
        sigx = 1.0e-10
    dtsft = (rhomn * ds * ds * ds / (sigx * 4.0 * (1.0 + cyl))) ** 0.5

    # --- RDTEXP ---
    rdtexp = 2.0 * abs(csq) ** 0.5 / ds
    if csq < 0.0:
        rdtexp = 1.0e10

    # --- BETA ---
    rcsq = 1.0 / (rhof * csq) if csq > 0.0 else 0.0
    ctos = delt * rdtexp
    comg = min(ctos * ctos, 1.0)
    omg1 = (omg - 1.0) * comg + 1.0

    for j in range(1, jm1):
        for i in range(1, im1):
            if BETA[i, j] < 0.0:
                continue

            dx_i = delx[i]
            dx_ip = delx[min(i + 1, ndelx - 1)]
            dx_im = delx[max(i - 1, 0)]
            dy_j = dely[j]
            dy_jp = dely[min(j + 1, ndely - 1)]
            dy_jm = dely[max(j - 1, 0)]

            rho_ij = rhofc + rhod * F[i, j]

            rhxr = rho_ij * dx_ip + (rhofc + rhod * F[i + 1, j]) * dx_i
            rhxl = (rhofc + rhod * F[i - 1, j]) * dx_i + rho_ij * dx_im
            rhyt = rho_ij * dy_jp + (rhofc + rhod * F[i, j + 1]) * dy_j
            rhyb = (rhofc + rhod * F[i, j - 1]) * dy_j + rho_ij * dy_jm

            xx = (delt * rdx[i] * (2.0 / rhxl + 2.0 / rhxr) +
                  delt * rdy[j] * (2.0 / rhyt + 2.0 / rhyb))

            rhor = rhof / rho_ij
            BETA[i, j] = omg1 / (xx * comg + rcsq * rhor / delt)

    return ds, dtvis, dtsft, rdtexp


def compute_beta(mesh, fields, state):
    """
    Compute BETA relaxation coefficients and stability limits.
    Port of Fortran SETUP lines 2641-2679.
    """
    ds, dtvis, dtsft, rdtexp = _compute_beta_kernel(
        fields.F, fields.BETA,
        mesh.delx, mesh.dely, mesh.rdx, mesh.rdy,
        mesh.im1, mesh.jm1, state.delt,
        state.rhof, state.rhofc, state.rhof - state.rhofc,
        state.csq, state.omg, state.vnu, state.vnuc, state.sigma,
        state.cyl)

    state.dtvis = dtvis
    state.dtsft = dtsft
    state.rdtexp = rdtexp


# =====================================================================
# DELTADJ — adaptive time stepping, port of Fortran DELTADJ
# =====================================================================

def deltadj(mesh, fields, state):
    """
    Adjust DELT based on CFL conditions, viscous/surface tension limits.
    Port of Fortran SUBROUTINE DELTADJ.
    """
    deltn = state.delt

    # --- Courant rollback ---
    if state.flgc >= 0.5:
        state.t -= state.delt
        state.icyle -= 1
        state.delt /= 2.0
        fields.P[:] = fields.PN
        fields.F[:] = fields.FN
        fields.U.fill(0.0)
        fields.V.fill(0.0)
        state.nflgc += 1

    # --- Auto time stepping ---
    if state.autot < 0.5 and state.fnoc < 0.5:
        # No auto-stepping, but still may need to recompute BETA
        if state.delt != deltn and state.nmat != 1:
            compute_beta(mesh, fields, state)
        return

    if state.fnoc > 0.5:
        state.delt /= 2.0

    # Find maximum velocity derivatives
    xi = mesh.xi        # cell centers
    yj = mesh.yj
    im1 = mesh.im1
    jm1 = mesh.jm1
    UN = fields.UN
    VN = fields.VN

    dumx = EM10
    dvmx = EM10

    nxi = len(xi)
    nyj = len(yj)

    for j in range(1, jm1):
        for i in range(1, im1):
            dxi = xi[min(i + 1, nxi - 1)] - xi[i] if i + 1 < nxi else xi[i] - xi[i - 1]
            dyj = yj[min(j + 1, nyj - 1)] - yj[j] if j + 1 < nyj else yj[j] - yj[j - 1]
            udm = abs(UN[i, j]) / dxi
            vdm = abs(VN[i, j]) / dyj
            dumx = max(dumx, udm)
            dvmx = max(dvmx, vdm)

    dtmp = 1.01
    if state.iter > 25:
        dtmp = 0.99
    delto = state.delt * dtmp
    con = 0.25
    delt_cfl_u = con / dumx
    delt_cfl_v = con / dvmx

    delt_new = min(delto, delt_cfl_u, delt_cfl_v, state.dtvis, state.dtsft)

    # Safety check
    delt_min = 1.0e-10
    if delt_new < delt_min:
        print(f"DELTADJ: time step too small ({delt_new:.3e}), "
              f"setting to {delt_min:.3e}")
        delt_new = delt_min

    state.delt = delt_new

    # Recompute BETA if DELT changed and NMAT != 1
    if state.delt != deltn and state.nmat != 1:
        compute_beta(mesh, fields, state)


# =====================================================================
# MPI-aware wrappers for compute_beta and deltadj
# =====================================================================

def compute_beta_mpi(mesh, fields, state, comm):
    """
    Compute BETA with global Allreduce for stability limits.
    Each rank computes local dtvis, dtsft, rdtexp, then we take
    the global MIN/MAX as appropriate.
    """
    ds, dtvis, dtsft, rdtexp = _compute_beta_kernel(
        fields.F, fields.BETA,
        mesh.delx, mesh.dely, mesh.rdx, mesh.rdy,
        mesh.im1, mesh.jm1, state.delt,
        state.rhof, state.rhofc, state.rhof - state.rhofc,
        state.csq, state.omg, state.vnu, state.vnuc, state.sigma,
        state.cyl)

    # Global reductions
    global_vals = np.array([dtvis, dtsft], dtype=float)
    comm.Allreduce(np.array([dtvis, dtsft], dtype=float),
                   global_vals, op=_mpi_min_op())
    state.dtvis = global_vals[0]
    state.dtsft = global_vals[1]

    # rdtexp depends on ds (min cell size globally)
    global_ds = np.array([ds], dtype=float)
    comm.Allreduce(np.array([ds], dtype=float), global_ds, op=_mpi_min_op())
    ds_global = global_ds[0]

    # Recompute rdtexp from global ds
    if state.csq >= 0.0:
        state.rdtexp = 2.0 * abs(state.csq) ** 0.5 / ds_global
    else:
        state.rdtexp = 1.0e10


def deltadj_mpi(mesh, fields, state, comm):
    """
    Adaptive time stepping with MPI reductions.
    """
    deltn = state.delt

    # --- Courant rollback (all ranks must agree) ---
    # First, Allreduce flgc so all ranks agree
    flgc_global = np.array([state.flgc], dtype=float)
    comm.Allreduce(np.array([state.flgc], dtype=float),
                   flgc_global, op=_mpi_max_op())
    state.flgc = flgc_global[0]

    if state.flgc >= 0.5:
        state.t -= state.delt
        state.icyle -= 1
        state.delt /= 2.0
        fields.P[:] = fields.PN
        fields.F[:] = fields.FN
        fields.U.fill(0.0)
        fields.V.fill(0.0)
        state.nflgc += 1

    # --- Auto time stepping ---
    if state.autot < 0.5 and state.fnoc < 0.5:
        if state.delt != deltn and state.nmat != 1:
            compute_beta_mpi(mesh, fields, state, comm)
        return

    if state.fnoc > 0.5:
        state.delt /= 2.0

    # Find maximum velocity derivatives (local)
    xi = mesh.xi
    yj = mesh.yj
    im1 = mesh.im1
    jm1 = mesh.jm1
    UN = fields.UN
    VN = fields.VN

    dumx = EM10
    dvmx = EM10

    nxi = len(xi)
    nyj = len(yj)

    for j in range(1, jm1):
        for i in range(1, im1):
            dxi = xi[min(i + 1, nxi - 1)] - xi[i] if i + 1 < nxi else xi[i] - xi[i - 1]
            dyj = yj[min(j + 1, nyj - 1)] - yj[j] if j + 1 < nyj else yj[j] - yj[j - 1]
            udm = abs(UN[i, j]) / dxi
            vdm = abs(VN[i, j]) / dyj
            dumx = max(dumx, udm)
            dvmx = max(dvmx, vdm)

    # Global max of velocity derivatives
    local_maxes = np.array([dumx, dvmx], dtype=float)
    global_maxes = np.empty(2, dtype=float)
    comm.Allreduce(local_maxes, global_maxes, op=_mpi_max_op())
    dumx = global_maxes[0]
    dvmx = global_maxes[1]

    dtmp = 1.01
    if state.iter > 25:
        dtmp = 0.99
    delto = state.delt * dtmp
    con = 0.25
    delt_cfl_u = con / dumx
    delt_cfl_v = con / dvmx

    delt_new = min(delto, delt_cfl_u, delt_cfl_v, state.dtvis, state.dtsft)

    delt_min = 1.0e-10
    if delt_new < delt_min:
        if comm.Get_rank() == 0:
            print(f"DELTADJ: time step too small ({delt_new:.3e}), "
                  f"setting to {delt_min:.3e}")
        delt_new = delt_min

    state.delt = delt_new

    if state.delt != deltn and state.nmat != 1:
        compute_beta_mpi(mesh, fields, state, comm)


def _mpi_min_op():
    from mpi4py import MPI
    return MPI.MIN


def _mpi_max_op():
    from mpi4py import MPI
    return MPI.MAX
