# bc.py
import numpy as np
from numba import njit
from .mesh import Mesh
from .fields import Fields
from .state import RunState
from .config import RunConfig


# ---------------------------------------------------------------------------
# Numba kernel for apply_boundary_conditions
# ---------------------------------------------------------------------------
@njit(cache=True)
def _apply_bc_kernel(U, V, P, F, PS, delx, dely, x, rx,
                     imax, jmax, iwl, iwr, iwt, iwb,
                     iter_, cyl, ishvel, shvel):
    """Numba-accelerated ghost-cell BC sweep."""
    fIM1 = imax - 2
    fIM2 = imax - 3
    fIMAX = imax - 1
    fJM1 = jmax - 2
    fJM2 = jmax - 3
    fJMAX = jmax - 1

    # Left/right boundaries
    for j in range(jmax):
        F[0, j] = F[1, j]
        F[fIMAX, j] = F[fIM1, j]
        P[0, j] = P[1, j]
        P[fIMAX, j] = P[fIM1, j]

        # Left wall
        if iwl == 1:
            U[0, j] = 0.0
            V[0, j] = V[1, j]
        elif iwl == 2:
            U[0, j] = 0.0
            V[0, j] = -V[1, j] * delx[0] / delx[1]
        elif iwl == 3 or iwl == 5:
            if iter_ <= 0:
                U[0, j] = U[1, j] * (x[1] * rx[0] * cyl + 1.0 - cyl)
                V[0, j] = V[1, j]
        elif iwl == 4:
            U[0, j] = U[fIM2, j]
            V[0, j] = V[fIM2, j]
            P[0, j] = P[fIM2, j]
            F[0, j] = F[fIM2, j]

        # Right wall
        if iwr == 1:
            U[fIM1, j] = 0.0
            V[fIMAX, j] = V[fIM1, j]
        elif iwr == 2:
            U[fIM1, j] = 0.0
            V[fIMAX, j] = -V[fIM1, j] * delx[fIMAX] / delx[fIM1]
        elif iwr == 3 or iwr == 5:
            if iter_ <= 0:
                U[fIM1, j] = U[fIM2, j] * (x[fIM2] * rx[fIM1] * cyl + 1.0 - cyl)
                V[fIMAX, j] = V[fIM1, j]
        elif iwr == 4:
            U[fIM1, j] = U[1, j]
            V[fIM1, j] = V[1, j]
            P[fIM1, j] = P[1, j]
            PS[fIM1, j] = PS[1, j]
            F[fIM1, j] = F[1, j]
            V[fIMAX, j] = V[2, j]
            F[fIMAX, j] = F[2, j]
            U[fIMAX, j] = U[2, j]

    # Bottom/top boundaries
    for i in range(imax):
        F[i, 0] = F[i, 1]
        F[i, fJMAX] = F[i, fJM1]
        P[i, 0] = P[i, 1]
        P[i, fJMAX] = P[i, fJM1]

        # Top wall
        if iwt == 1:
            V[i, fJM1] = 0.0
            U[i, fJMAX] = U[i, fJM1]
        elif iwt == 2:
            V[i, fJM1] = 0.0
            U[i, fJMAX] = -U[i, fJM1] * dely[fJMAX] / dely[fJM1]
        elif iwt == 3 or iwt == 5:
            if iter_ <= 0:
                V[i, fJM1] = V[i, fJM2]
                U[i, fJMAX] = U[i, fJM1]
        elif iwt == 4:
            V[i, fJM1] = V[i, 1]
            U[i, fJM1] = U[i, 1]
            P[i, fJM1] = P[i, 1]
            PS[i, fJM1] = PS[i, 1]
            F[i, fJM1] = F[i, 1]
            U[i, fJMAX] = U[i, 2]
            F[i, fJMAX] = F[i, 2]
            V[i, fJMAX] = V[i, 2]

        # Bottom wall
        if iwb == 1:
            V[i, 0] = 0.0
            U[i, 0] = U[i, 1]
        elif iwb == 2:
            V[i, 0] = 0.0
            U[i, 0] = -U[i, 1]
        elif iwb == 3 or iwb == 5:
            if iter_ <= 0:
                V[i, 0] = V[i, 1]
                U[i, 0] = U[i, 1]
        elif iwb == 4:
            V[i, 0] = V[i, fJM2]
            U[i, 0] = U[i, fJM2]
            F[i, 0] = F[i, fJM2]

    # Shear velocity BC
    if ishvel:
        for i in range(imax):
            U[i, 0] = shvel
            U[i, 1] = shvel
            U[i, fJMAX] = -shvel
            U[i, fJM1] = -shvel


def apply_boundary_conditions(mesh: Mesh, fields: Fields, state: RunState,
                              iwl: int, iwr: int, iwt: int, iwb: int) -> None:
    """
    Python version of SUBROUTINE BC.
    Fills ghost cells of U, V, P, F based on wall types and periodicity.

    Index mapping (Fortran 1-based -> Python 0-based):
      Fortran I=1    -> Python 0
      Fortran I=2    -> Python 1
      Fortran I=IM1  -> Python imax-2  (last interior cell)
      Fortran I=IM2  -> Python imax-3
      Fortran I=IMAX -> Python imax-1  (right ghost cell)
    Note: mesh.im1 = imax-1 = Python index of Fortran IMAX, NOT IM1!
    """
    _apply_bc_kernel(
        fields.U, fields.V, fields.P, fields.F, fields.PS,
        mesh.delx, mesh.dely, mesh.x, mesh.rx,
        mesh.imax, mesh.jmax, iwl, iwr, iwt, iwb,
        state.iter, state.cyl,
        1 if state.ishvel else 0, state.shvel)


# ---------------------------------------------------------------------------
# Numba kernel for free_surface_bc
# ---------------------------------------------------------------------------
EMF = 1.0e-6
EMF1 = 1.0 - EMF


@njit(cache=True)
def _free_surface_bc_kernel(U, V, P, F, BETA,
                            rdx, rdy, rxi, delx, dely,
                            i_lo, i_hi, j_lo, j_hi,
                            cyl, nmat, flg, iter_):
    """
    Numba-accelerated free-surface and obstacle BC.

    Part 1: Obstacle cell averaging (BETA<0).
    Part 2: NMAT=1 surface cell velocity extrapolation (NFSB bitmask).
    Part 3: Propagate velocities into empty neighbours.
    """
    emf = 1.0e-6
    emf1 = 1.0 - emf

    for i in range(i_lo, i_hi):
        xrp = rdx[i] + rxi[i] / 2.0
        rxrp = 1.0 / xrp if xrp != 0.0 else 0.0
        xrm = rdx[i] - rxi[i] / 2.0
        rxrm = 1.0 / xrm if xrm > 0.0 else 0.0

        for j in range(j_lo, j_hi):
            # Part 1: obstacle cell averaging
            if BETA[i, j] < 0.0:
                bmr = 1.0 if BETA[i + 1, j] > 0.0 else 0.0
                bmt = 1.0 if BETA[i, j + 1] > 0.0 else 0.0
                bml = 1.0 if BETA[i - 1, j] > 0.0 else 0.0
                bmb = 1.0 if BETA[i, j - 1] > 0.0 else 0.0
                bmtot = bmr + bmt + bml + bmb
                F[i, j] = 0.0
                P[i, j] = 0.0
                if bmtot > 0.0:
                    F[i, j] = (bmr * F[i + 1, j] + bmt * F[i, j + 1] +
                               bml * F[i - 1, j] + bmb * F[i, j - 1]) / bmtot
                    P[i, j] = (bmr * P[i + 1, j] + bmt * P[i, j + 1] +
                               bml * P[i - 1, j] + bmb * P[i, j - 1]) / bmtot
                continue

            # Part 2: NMAT=1 free-surface velocity extrapolation
            if nmat == 2:
                continue
            if F[i, j] < emf or F[i, j] > emf1:
                continue

            nfsb = 0
            if F[i + 1, j] < emf:
                nfsb += 1
            if F[i, j + 1] < emf:
                nfsb += 2
            if F[i - 1, j] < emf:
                nfsb += 4
            if F[i, j - 1] < emf:
                nfsb += 8
            if nfsb == 0:
                continue

            if nfsb == 1:
                U[i, j] = (U[i - 1, j] - delx[i] * rdy[j] * (V[i, j] - V[i, j - 1])) * (1.0 - cyl) + \
                           cyl * (U[i - 1, j] * xrm * rxrp - rdy[j] * rxrp * (V[i, j] - V[i, j - 1]))
            elif nfsb == 2:
                V[i, j] = (V[i, j - 1] - dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                           cyl * (V[i, j - 1] - dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 3:
                U[i, j] = U[i - 1, j]
                V[i, j] = (V[i, j - 1] - dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                           cyl * (V[i, j - 1] - dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 4:
                U[i - 1, j] = (U[i, j] + delx[i] * rdy[j] * (V[i, j] - V[i, j - 1])) * (1.0 - cyl) + \
                               cyl * (U[i, j] * xrp * rxrm + rdy[j] * rxrm * (V[i, j] - V[i, j - 1]))
            elif nfsb == 5:
                U[i - 1, j] = U[i - 1, j - 1]
                U[i, j] = (U[i - 1, j] - delx[i] * rdy[j] * (V[i, j] - V[i, j - 1])) * (1.0 - cyl) + \
                           cyl * (U[i - 1, j] * xrm * rxrp - rdy[j] * rxrp * (V[i, j] - V[i, j - 1]))
            elif nfsb == 6:
                U[i - 1, j] = U[i, j] * (1.0 - cyl) + cyl * U[i, j]
                V[i, j] = (V[i, j - 1] - dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                           cyl * (V[i, j - 1] - dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 7:
                U[i - 1, j] = U[i - 1, j - 1]
                U[i, j] = U[i, j - 1]
                V[i, j] = (V[i, j - 1] - dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                           cyl * (V[i, j - 1] - dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 8:
                V[i, j - 1] = (V[i, j] + dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                               cyl * (V[i, j] + dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 9:
                U[i, j] = U[i - 1, j] * (1.0 - cyl) + cyl * U[i - 1, j]
                V[i, j - 1] = (V[i, j] + dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                               cyl * (V[i, j] + dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 10:
                V[i, j] = V[i - 1, j]
                V[i, j - 1] = (V[i, j] + dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                               cyl * (V[i, j] + dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 11:
                V[i, j] = V[i - 1, j]
                V[i, j - 1] = V[i - 1, j - 1]
                U[i, j] = (U[i - 1, j] - delx[i] * rdy[j] * (V[i, j] - V[i, j - 1])) * (1.0 - cyl) + \
                           cyl * (U[i - 1, j] * xrm * rxrp - rdy[j] * rxrp * (V[i, j] - V[i, j - 1]))
            elif nfsb == 12:
                U[i - 1, j] = U[i, j] * (1.0 - cyl) + cyl * U[i, j]
                V[i, j - 1] = (V[i, j] + dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                               cyl * (V[i, j] + dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 13:
                U[i, j] = U[i, j + 1]
                U[i - 1, j] = U[i - 1, j + 1]
                V[i, j - 1] = (V[i, j] + dely[j] * rdx[i] * (U[i, j] - U[i - 1, j])) * (1.0 - cyl) + \
                               cyl * (V[i, j] + dely[j] * (xrp * U[i, j] - xrm * U[i - 1, j]))
            elif nfsb == 14:
                V[i, j] = V[i + 1, j]
                V[i, j - 1] = V[i + 1, j - 1]
                U[i - 1, j] = (U[i, j] + delx[i] * rdy[j] * (V[i, j] - V[i, j - 1])) * (1.0 - cyl) + \
                               cyl * (U[i, j] * xrp * rxrm + rdy[j] * rxrm * (V[i, j] - V[i, j - 1]))
            elif nfsb == 15:
                U[i, j] = U[i - 1, j] * (1.0 - cyl) + cyl * U[i - 1, j] * xrm * rxrp
                V[i, j - 1] = V[i, j]
                V[i, j + 1] = V[i, j]

            # Part 3: set velocities in empty cells adjacent to surface
            if flg > 0.5 and iter_ > 0:
                continue

            if F[i + 1, j] < emf:
                if F[i + 1, j + 1] < emf:
                    V[i + 1, j] = V[i, j]
                if F[i + 1, j - 1] < emf:
                    V[i + 1, j - 1] = V[i, j - 1]
            if F[i, j + 1] < emf:
                if F[i + 1, j + 1] < emf:
                    U[i, j + 1] = U[i, j]
                if F[i - 1, j + 1] < emf:
                    U[i - 1, j + 1] = U[i - 1, j]
            if F[i - 1, j] < emf:
                if F[i - 1, j + 1] < emf:
                    V[i - 1, j] = V[i, j]
                if F[i - 1, j - 1] < emf:
                    V[i - 1, j - 1] = V[i, j - 1]
            if F[i, j - 1] < emf:
                if F[i + 1, j - 1] < emf:
                    U[i, j - 1] = U[i, j]
                if F[i - 1, j - 1] < emf:
                    U[i - 1, j - 1] = U[i - 1, j]


def free_surface_bc(mesh: Mesh, fields: Fields, state: RunState,
                    cfg: RunConfig) -> None:
    """
    Port of Fortran BC lines 369-466: free-surface and obstacle boundary conditions.

    Part 1: Obstacle cell averaging (BETA<0): average F,P from fluid neighbors.
    Part 2: NMAT=1 surface cell velocity extrapolation using NFSB bitmask.
    Part 3: Propagate velocities into empty neighbors (label 410).
    """
    _free_surface_bc_kernel(
        fields.U, fields.V, fields.P, fields.F, fields.BETA,
        mesh.rdx, mesh.rdy, mesh.rxi, mesh.delx, mesh.dely,
        1, mesh.imax - 2 + 1, 1, mesh.jmax - 2 + 1,
        state.cyl, state.nmat, state.flg, state.iter)
