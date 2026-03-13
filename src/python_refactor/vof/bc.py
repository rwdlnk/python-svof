# bc.py
import numpy as np
from .mesh import Mesh
from .fields import Fields
from .state import RunState


def apply_boundary_conditions(mesh: Mesh, fields: Fields, state: RunState,
                              iwl: int, iwr: int, iwt: int, iwb: int) -> None:
    """
    Python version of SUBROUTINE BC.[file:40]
    Fills ghost cells of U, V, P, F based on wall types and periodicity.
    """
    U, V, P, F = fields.U, fields.V, fields.P, fields.F
    imax, jmax = mesh.imax, mesh.jmax
    im1, jm1 = mesh.im1, mesh.jm1
    im2, jm2 = mesh.im2, mesh.jm2
    delx, dely = mesh.delx, mesh.dely
    x = mesh.x
    rx = mesh.rx
    rdx = mesh.rdx
    rdy = mesh.rdy

    # Left/right boundaries: loop over J
    for j in range(jmax):
        # copy F and P at domain ends
        F[0, j] = F[1, j]
        F[imax - 1, j] = F[im1, j]
        P[0, j] = P[1, j]
        P[imax - 1, j] = P[im1, j]

        # Left wall I=0 (Fortran I=1)
        if iwl == 1:
            # rigid free-slip
            U[0, j] = 0.0
            V[0, j] = V[1, j]
        elif iwl == 2:
            # rigid no-slip
            U[0, j] = 0.0
            V[0, j] = -V[1, j] * delx[0] / delx[1]
        elif iwl in (3, 5):
            # continuance or specified pressure; use interior values only for iter<=0
            if state.iter <= 0:
                U[0, j] = U[1, j] * (x[1] * rx[0] * state.epsi + 1.0 - state.epsi)
                V[0, j] = V[1, j]
        elif iwl == 4:
            # periodic in x: copy from right side interior
            U[0, j] = U[im2, j]
            V[0, j] = V[im2, j]
            F[0, j] = F[im2, j]

        # Right wall I=imax-1 (Fortran I=IMAX)
        if iwr == 1:
            U[im1, j] = 0.0
            V[imax - 1, j] = V[im1, j]
        elif iwr == 2:
            U[im1, j] = 0.0
            V[imax - 1, j] = -V[im1, j] * delx[imax - 1] / delx[im1]
        elif iwr in (3, 5):
            if state.iter <= 0:
                U[im1, j] = U[im2, j] * (x[im2] * rx[im1] * state.epsi + 1.0 - state.epsi)
                V[imax - 1, j] = V[im1, j]
        elif iwr == 4:
            U[im1, j] = U[1, j]
            V[im1, j] = V[1, j]
            P[im1, j] = P[1, j]
            F[im1, j] = F[1, j]
            V[imax - 1, j] = V[2, j]
            F[imax - 1, j] = F[2, j]

    # Bottom/top boundaries: loop over I
    for i in range(imax):
        F[i, 0] = F[i, 1]
        F[i, jmax - 1] = F[i, jm1]
        P[i, 0] = P[i, 1]
        P[i, jmax - 1] = P[i, jm1]

        # Top wall (Fortran J=JMAX)
        if iwt == 1:
            V[i, jm1] = 0.0
            U[i, jmax - 1] = U[i, jm1]
        elif iwt == 2:
            V[i, jm1] = 0.0
            U[i, jmax - 1] = -U[i, jm1] * dely[jmax - 2] / dely[jm1-1]
        elif iwt in (3, 5):
            if state.iter <= 0:
                V[i, jm1] = V[i, jm2]
                U[i, jmax - 1] = U[i, jm1]
        elif iwt == 4:
            V[i, jm1] = V[i, 1]
            U[i, jm1] = U[i, 1]
            P[i, jm1] = P[i, 1]
            F[i, jm1] = F[i, 1]
            U[i, jmax - 1] = U[i, 2]
            F[i, jmax - 1] = F[i, 2]

        # Bottom wall (Fortran J=1)
        if iwb == 1:
            V[i, 0] = 0.0
            U[i, 0] = U[i, 1]
        elif iwb == 2:
            V[i, 0] = 0.0
            U[i, 0] = -U[i, 1]
        elif iwb in (3, 5):
            if state.iter <= 0:
                V[i, 0] = V[i, 1]
                U[i, 0] = U[i, 1]
        elif iwb == 4:
            V[i, 0] = V[i, jm2]
            U[i, 0] = U[i, jm2]
            F[i, 0] = F[i, jm2]

    # Special shear-velocity BC (ISHVEL)
    if state.ishvel:
        for i in range(imax):
            U[i, 0] = state.shvel
            U[i, 1] = state.shvel
            U[i, jmax - 1] = -state.shvel
            U[i, jm1] = -state.shvel

