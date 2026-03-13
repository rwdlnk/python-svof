# initcond.py
import numpy as np
from .mesh import Mesh
from .fields import Fields
from .config import RunConfig
from .state import RunState


def setup_initial_state(mesh: Mesh, fields: Fields,
                        cfg: RunConfig, state: RunState) -> None:
    """
    Rough Python analogue of SETUP.[file:40]
    For now:
      - zero U,V,P
      - initialize F as a horizontal interface at Y = cfg.h_dist (if NMAT=2)
      - zero PS, NF, BETA.
    """
    U, V, P, F = fields.U, fields.V, fields.P, fields.F
    PS, NF, BETA = fields.PS, fields.NF, fields.BETA

    U.fill(0.0)
    V.fill(0.0)
    P.fill(0.0)
    PS.fill(0.0)
    NF.fill(0)
    BETA.fill(1.0)  # all fluid, no obstacles yet

    if cfg.nmat == 2:
        # Simple two-fluid setup: F=1 below interface, 0 above
        yj = mesh.yj
        h = cfg.h_dist
        for j in range(mesh.jmax-1):
            if yj[j] <= h:
                F[:, j] = 1.0
            else:
                F[:, j] = 0.0
    else:
        F.fill(1.0)

    # Hydrostatic pressure: integrate from top down
    # Mirrors Fortran: P(I,JMAX)=0, P(I,JJ)=P(I,JJ+1)-GY*RHOYA
    gy = cfg.gy
    rhof = cfg.rhof
    rhofc = cfg.rhofc
    rhod = rhof - rhofc
    nmat = cfg.nmat
    dely = mesh.dely
    im1, jm1 = mesh.im1, mesh.jm1
    jmax = mesh.jmax

    # P=0 at top boundary (Python index jmax-1 = Fortran JMAX)
    P[:, jmax - 1] = 0.0
    # Integrate downward: Fortran J=2..JM1 with JJ=JM1-J+2 (JM1 down to 2)
    # In 0-based: jj goes from jm1-1 down to 1 (jj+1 from jm1 down to 2)
    ndely = len(dely)
    for i in range(1, im1):
        for jj in range(jm1 - 1, 0, -1):
            djj  = dely[jj]
            djj1 = dely[min(jj + 1, ndely - 1)]   # clamp for top ghost
            if nmat == 2:
                rhoya = ((rhofc + rhod * F[i, jj]) * djj / 2.0 +
                         (rhofc + rhod * F[i, jj + 1]) * djj1 / 2.0)
            else:
                rhoya = (min(F[i, jj + 1], 0.5) * djj1 +
                         max(0.0, F[i, jj] - 0.5) * djj) * rhof
            P[i, jj] = P[i, jj + 1] - gy * rhoya

    # Copy hydrostatic P to ghost cells
    P[0, :] = P[1, :]
    P[im1, :] = P[im1 - 1, :]
    P[:, 0] = P[:, 1]

    # Apply initial disturbance if requested
    apply_disturbance(mesh, fields, cfg, state)

def apply_disturbance(mesh: Mesh, fields: Fields,
                      cfg: RunConfig, state: RunState) -> None:
    """
    Dispatcher mirroring Fortran DISTURB (NDIS, IDIS, VORIG, R,H,ILOW,IHIGH).[file:40][file:57]
    """
    if cfg.ndis == 0:
        return
    elif cfg.ndis == 1:
        _planar_disturbance(mesh, fields, cfg)
    elif cfg.ndis == 2:
        # TODO: cylindrical disturbance (CYLIND)
        _planar_disturbance(mesh, fields, cfg)
    elif cfg.ndis == 3:
        _noise_disturbance(mesh, fields, cfg)
    elif cfg.ndis == 4:
        # TODO: PLANA2
        _planar_disturbance(mesh, fields, cfg)
    elif cfg.ndis in (5, 6):
        # TODO: RTINIT / DALZIEL_PERTURB
        _planar_disturbance(mesh, fields, cfg)
    else:
        # Unknown NDIS: do nothing but warn
        print(f"Warning: NDIS={cfg.ndis} is invalid; no disturbance applied.")

def _noise_disturbance(mesh: Mesh, fields: Fields, cfg: RunConfig) -> None:
    """
    Python version of NOISE: white-noise initial velocities scaled by VORIG.[file:40]
    """
    U, V = fields.U, fields.V
    F = fields.F
    jmax = mesh.jmax
    imax = mesh.imax

    rng = np.random.default_rng()
    for i in range(imax):
        for j in range(jmax):
            # Skip solid/void cells later once BETA/obstacles are wired in
            uamp = (rng.random() - 0.5) * cfg.vorig
            vamp = (rng.random() - 0.5) * cfg.vorig
            if abs(uamp) >= 1.0e-10:
                U[i, j] += uamp
            if abs(vamp) >= 1.0e-10:
                V[i, j] += vamp


def _planar_disturbance(mesh: Mesh, fields: Fields, cfg: RunConfig) -> None:
    """
    Port of Fortran SUBROUTINE PLANAR.
    Rayleigh-Taylor disturbance per Youngs, Physica 12D (1984) 32-44.

    Accumulates cosine modes k = ILOW..IHIGH with random coefficients
    in [-1, +1], then scales the total by VORIG so that VORIG controls
    the maximum disturbance amplitude.

    V is uniform in y (set for all rows j = 0..JM1-1).
    """
    U, V = fields.U, fields.V
    xi = mesh.xi        # cell centers, length imax-1
    imax = mesh.imax
    jm1 = mesh.jm1

    r = cfg.r_dist      # x-width of disturbance region
    ilow = cfg.ilow
    ihigh = cfg.ihigh
    vorig = cfg.vorig
    pi = np.pi

    # Zero U and V (Fortran does this at the top of PLANAR)
    U.fill(0.0)
    V.fill(0.0)

    # Use a fixed seed for reproducibility across runs
    rng = np.random.default_rng(seed=42)

    # Accumulate modes with random coefficients in [-1, +1]
    nxi = len(xi)       # imax - 1
    for k in range(ilow, ihigh + 1):
        coeff = 2.0 * rng.random() - 1.0       # uniform in [-1, +1]
        pos = 2.0 * k * pi * xi[:nxi] / r
        vsum = coeff * np.cos(pos)              # shape (nxi,)
        V[:nxi, 0:jm1] += vsum[:, np.newaxis]

    # Scale total by VORIG
    V[:nxi, 0:jm1] *= vorig

