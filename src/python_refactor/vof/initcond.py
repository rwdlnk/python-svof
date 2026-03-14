# initcond.py
import numpy as np
from .mesh import Mesh
from .fields import Fields
from .config import RunConfig
from .state import RunState

# Coefficients from gfortran RAND(0) with default seed, transformed as 2*RAND(0)-1.
# These match the Fortran reference exactly for modes k=1..10 (ILOW..IHIGH).
_PLANAR_COEFFS = [
    2.0 * 0.00000762939453125 - 1.0,    # -0.99999237060546875
    2.0 * 0.13153767585754395 - 1.0,    # -0.73692464828491211
    2.0 * 0.75560522079467773 - 1.0,    #  0.51121044158935547
    2.0 * 0.45865011215209961 - 1.0,    # -0.08269977569580078
    2.0 * 0.53276705741882324 - 1.0,    #  0.06553411483764648
    2.0 * 0.21895909309387207 - 1.0,    # -0.56208181381225586
    2.0 * 0.04704451560974121 - 1.0,    # -0.90591096878051758
    2.0 * 0.67886447906494141 - 1.0,    #  0.35772895812988281
    2.0 * 0.67929625511169434 - 1.0,    #  0.35859251022338867
    2.0 * 0.93469285964965820 - 1.0,    #  0.86938571929931641
]


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
    BETA.fill(1.0)  # all fluid initially

    # Set BETA = -1 in obstacle cells (Fortran: SETUP lines 2681-2694)
    if cfg.nobs > 0:
        # Convert 1-based Fortran indices to 0-based Python
        iomin0 = cfg.iomin - 1
        iomax0 = cfg.iomax      # exclusive upper bound
        jomin0 = cfg.jomin - 1
        jomax0 = cfg.jomax      # exclusive upper bound
        BETA[iomin0:iomax0, jomin0:jomax0] = -1.0

    # F field initialization uses IFL/IFR/JFB/JFT for both NMAT=1 and NMAT=2.
    # Fortran: F=0 everywhere, then F=1 for I=IFL..IFR, J=JFB..JFT.
    # Convert to 0-based: i=ifl-1..ifr-1, j=jfb-1..jft-1.
    F.fill(0.0)
    ifl0 = cfg.ifl - 1
    ifr0 = cfg.ifr       # exclusive upper bound in 0-based
    jfb0 = cfg.jfb - 1
    jft0 = cfg.jft        # exclusive upper bound in 0-based
    F[ifl0:ifr0, jfb0:jft0] = 1.0

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
    Dispatcher mirroring Fortran DISTURB (NDIS, IDIS, VORIG, R,H,ILOW,IHIGH).
    """
    if cfg.ndis == 0:
        return
    elif cfg.ndis == 1:
        _planar_disturbance(mesh, fields, cfg)
    elif cfg.ndis == 2:
        print("Warning: NDIS=2 (CYLIND) not implemented; requires Bessel functions.")
        print("         Falling back to PLANAR disturbance.")
        _planar_disturbance(mesh, fields, cfg)
    elif cfg.ndis == 3:
        _noise_disturbance(mesh, fields, cfg)
    elif cfg.ndis == 4:
        _plana2_disturbance(mesh, fields, cfg)
    elif cfg.ndis == 5:
        _rtinit_disturbance(mesh, fields, cfg)
    elif cfg.ndis == 6:
        _dalziel_perturb(mesh, fields, cfg)
    else:
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

    Accumulates cosine modes k = ILOW..IHIGH with fixed coefficients
    extracted from gfortran RAND(0) default seed, then scales the total
    by VORIG so that VORIG controls the maximum disturbance amplitude.

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

    # Accumulate modes with fixed coefficients matching Fortran RAND(0)
    nxi = len(xi)       # imax - 1
    for idx, k in enumerate(range(ilow, ihigh + 1)):
        coeff = _PLANAR_COEFFS[idx]
        pos = 2.0 * k * pi * xi[:nxi] / r
        vsum = coeff * np.cos(pos)              # shape (nxi,)
        V[:nxi, 0:jm1] += vsum[:, np.newaxis]

    # Scale total by VORIG
    V[:nxi, 0:jm1] *= vorig


def _plana2_disturbance(mesh: Mesh, fields: Fields, cfg: RunConfig) -> None:
    """
    Port of Fortran SUBROUTINE PLANA2 (lines 1708-1769).
    Rayleigh-Taylor disturbance per Youngs (1984) for planar geometry.
    Uses AK = I*PI/H (not 2*K*PI/R like PLANAR).
    """
    U, V = fields.U, fields.V
    BETA = fields.BETA
    x = mesh.x       # cell-face x coordinates
    y = mesh.y       # cell-face y coordinates
    imax = mesh.imax
    jmax = mesh.jmax

    r = cfg.r_dist
    h = cfg.h_dist
    ilow = cfg.ilow
    ihigh = cfg.ihigh
    pi = np.pi

    U.fill(0.0)
    V.fill(0.0)

    g = np.sqrt(cfg.gx ** 2 + cfg.gy ** 2)
    s = cfg.rhof / cfg.rhofc if cfg.rhofc != 0.0 else 1.0
    ahigh = float(ihigh - ilow + 1)

    # Use gfortran RAND(0) sequence for reproducibility
    rng = _FortranRand()

    for ii in range(imax):
        xr = x[ii]
        for jj in range(jmax):
            if BETA[ii, jj] < 0.0:
                continue

            usum = 0.0
            vsum = 0.0
            z = y[jj] - h
            if z == 0.0:
                z = 1.0e-8

            for i in range(ilow, ihigh + 1):
                randum = rng.next()
                amp = ((randum - 0.5) * 2.0e-3) / ahigh
                ar = r * amp
                ak = float(i) * pi / h
                usum += ar * np.cos(ak * xr) * np.exp(-ak * z) * z / abs(z)
                vsum += ar * ak * np.sin(ak * xr) * np.exp(-ak * z)

            vsum = -vsum
            if abs(usum) >= 1.0e-10:
                U[ii, jj] = usum
            if abs(vsum) >= 1.0e-10:
                V[ii, jj] = -vsum


class _FortranRand:
    """Mimics gfortran RAND(0) using default seed sequence."""
    # Pre-computed first 100 values from gfortran RAND(0)
    # For modes beyond the pre-computed set, we use a simple LCG
    _PRECOMPUTED = [
        0.00000762939453125, 0.13153767585754395, 0.75560522079467773,
        0.45865011215209961, 0.53276705741882324, 0.21895909309387207,
        0.04704451560974121, 0.67886447906494141, 0.67929625511169434,
        0.93469285964965820,
    ]

    def __init__(self):
        self._idx = 0
        # LCG state for values beyond precomputed
        self._lcg_state = 123456789

    def next(self):
        if self._idx < len(self._PRECOMPUTED):
            val = self._PRECOMPUTED[self._idx]
            self._idx += 1
            return val
        # Fallback LCG
        self._lcg_state = (self._lcg_state * 1103515245 + 12345) % (2 ** 31)
        self._idx += 1
        return self._lcg_state / (2 ** 31)


def _rtinit_disturbance(mesh: Mesh, fields: Fields, cfg: RunConfig) -> None:
    """
    Port of Fortran SUBROUTINE RTINIT (lines 2066-2434).
    Linear stability theory initialization for Rayleigh-Taylor.
    Supports periodic/non-periodic lateral BCs and free-slip/no-slip vertical BCs.
    """
    U, V = fields.U, fields.V
    imax = mesh.imax
    jmax = mesh.jmax
    xi = mesh.xi
    yj = mesh.yj
    delx = mesh.delx

    hb = cfg.h_dist
    l = cfg.r_dist
    nummodes = cfg.ihigh - cfg.ilow

    wl, wr, wb, wt = cfg.iwl, cfg.iwr, cfg.iwb, cfg.iwt
    vorig = cfg.vorig

    pi = np.pi
    twopi = 2.0 * pi
    alpha_pow = 1.5
    a0 = 0.01

    # Domain height above interface
    ht = mesh.y[mesh.ibar + 1] - mesh.y[1] - hb if mesh.ibar + 1 < len(mesh.y) else mesh.y[-1] - mesh.y[1] - hb

    drho = cfg.rhof - cfg.rhofc
    dx = delx[1]  # Fortran DELX(2)
    g = np.sqrt(cfg.gx ** 2 + cfg.gy ** 2)

    # Effective viscosity
    nueff = (cfg.rhofc * cfg.vnuc + cfg.rhof * cfg.vnu) / (cfg.rhofc + cfg.rhof) \
        if (cfg.rhofc + cfg.rhof) > 0 else 0.0
    inviscid = nueff <= 1.0e-6

    lateral_periodic = (wl == 4 and wr == 4)
    vertical_noslip = (wb == 2 and wt == 2)

    print(f"RTINIT: BCs: WL={wl} WR={wr} WB={wb} WT={wt}")
    print(f"  Lateral periodic: {lateral_periodic}")
    print(f"  Vertical no-slip: {vertical_noslip}")

    # Calculate wavenumber parameters
    if inviscid:
        if cfg.sigma > 0.0:
            kc = np.sqrt(g * abs(drho) / cfg.sigma)
            kmax = np.sqrt(g * abs(drho) / (3.0 * cfg.sigma))
        else:
            kc = pi / dx
            kmax = pi / (10.0 * dx)
    else:
        if cfg.sigma > 0.0:
            kc = np.sqrt(g * abs(drho) / cfg.sigma)
            kc = np.sqrt(g * abs(drho) / (cfg.sigma +
                         2.0 * nueff ** 2 * (cfg.rhofc + cfg.rhof) ** 2 / kc))
            kmax = np.sqrt(g * abs(drho) / (3.0 * cfg.sigma +
                           4.0 * nueff ** 2 * (cfg.rhofc + cfg.rhof) ** 2 / cfg.sigma))
        else:
            knu = np.sqrt(g * abs(drho) / (2.0 * nueff * (cfg.rhofc + cfg.rhof)))
            kc = knu
            kmax = knu / np.sqrt(3.0)

    # Number of modes
    if lateral_periodic:
        nmax = int(kc * l / twopi)
    else:
        nmax = int(kc * l / pi)
    nmax = min(nmax, nummodes)

    nx = mesh.ibar + 2
    ny = mesh.jbar + 2

    U.fill(0.0)
    V.fill(0.0)

    # Add contributions from each mode
    for n in range(1, nmax + 1):
        if lateral_periodic:
            kx = n * twopi / l
        else:
            kx = n * pi / l

        k2 = kx * kx
        an = a0 * (1.0 / n) ** alpha_pow

        # Pre-compute hyperbolic functions for no-slip
        if vertical_noslip:
            sinh_kb = np.sinh(kx * hb)
            sinh_kt = np.sinh(kx * ht)
            cosh_kb = np.cosh(kx * hb)
            cosh_kt = np.cosh(kx * ht)

        # Loop through interior (Fortran 2..NY-1 → Python 1..ny-2)
        for j in range(1, ny - 1):
            z = yj[j] - hb

            for i in range(1, nx - 1):
                x_coord = xi[i]

                # Horizontal structure
                if lateral_periodic:
                    phase_x_w = np.cos(kx * x_coord)
                    phase_x_u = np.cos(kx * x_coord)
                elif wl == 1 and wr == 1:
                    phase_x_w = np.cos(kx * x_coord)
                    phase_x_u = np.sin(kx * x_coord)
                else:
                    phase_x_w = np.cos(kx * x_coord)
                    phase_x_u = np.cos(kx * x_coord)

                if z < 0.0:
                    # Lower fluid
                    if vertical_noslip:
                        phi1 = np.sinh(kx * (z + hb)) / sinh_kb if abs(sinh_kb) > 1e-30 else 0.0
                        dphi1 = kx * np.cosh(kx * (z + hb)) / sinh_kb if abs(sinh_kb) > 1e-30 else 0.0
                    elif wb == 1 and wt == 1:
                        exp_term = np.exp(-2.0 * kx * hb)
                        denom = 1.0 - exp_term
                        if abs(denom) > 1e-30:
                            phi1 = (np.exp(-kx * abs(z)) - np.exp(kx * abs(z) - 2.0 * kx * hb)) / denom
                            dphi1 = kx * (np.exp(-kx * abs(z)) + np.exp(kx * abs(z) - 2.0 * kx * hb)) / denom
                        else:
                            phi1 = 0.0
                            dphi1 = 0.0
                    else:
                        sinh_khb = np.sinh(kx * hb)
                        phi1 = np.sinh(kx * (z + hb)) / sinh_khb if abs(sinh_khb) > 1e-30 else 0.0
                        dphi1 = kx * np.cosh(kx * (z + hb)) / sinh_khb if abs(sinh_khb) > 1e-30 else 0.0

                    U[i, j] += an * (-kx / np.sqrt(k2) * dphi1) * phase_x_u
                    V[i, j] += an * phi1 * phase_x_w
                else:
                    # Upper fluid
                    if vertical_noslip:
                        phi2 = np.sinh(kx * (ht - z)) / sinh_kt if abs(sinh_kt) > 1e-30 else 0.0
                        dphi2 = -kx * np.cosh(kx * (ht - z)) / sinh_kt if abs(sinh_kt) > 1e-30 else 0.0
                    elif wb == 1 and wt == 1:
                        exp_term = np.exp(-2.0 * kx * ht)
                        denom = 1.0 - exp_term
                        if abs(denom) > 1e-30:
                            phi2 = (np.exp(-kx * z) - np.exp(kx * z - 2.0 * kx * ht)) / denom
                            dphi2 = -kx * (np.exp(-kx * z) + np.exp(kx * z - 2.0 * kx * ht)) / denom
                        else:
                            phi2 = 0.0
                            dphi2 = 0.0
                    else:
                        sinh_kht = np.sinh(kx * ht)
                        phi2 = np.sinh(kx * (ht - z)) / sinh_kht if abs(sinh_kht) > 1e-30 else 0.0
                        dphi2 = -kx * np.cosh(kx * (ht - z)) / sinh_kht if abs(sinh_kht) > 1e-30 else 0.0

                    U[i, j] += an * (-kx / np.sqrt(k2) * dphi2) * phase_x_u
                    V[i, j] += an * phi2 * phase_x_w

    # Cap velocities at VORIG
    if vorig > 0.0:
        for j in range(1, ny - 1):
            for i in range(1, nx - 1):
                U[i, j] = max(-vorig, min(vorig, U[i, j]))
                V[i, j] = max(-vorig, min(vorig, V[i, j]))
        print(f"RTINIT: Velocities capped at VORIG = {vorig}")

    # Apply ghost cell BCs
    for j in range(ny):
        if wl == 1:
            U[0, j] = 0.0
            V[0, j] = V[1, j]
        elif wl == 2:
            U[0, j] = -U[1, j]
            V[0, j] = -V[1, j]
        elif wl in (3, 5):
            U[0, j] = U[1, j]
            V[0, j] = V[1, j]
        elif wl == 4:
            U[0, j] = U[nx - 2, j]
            V[0, j] = V[nx - 2, j]

        if wr == 1:
            U[nx - 1, j] = 0.0
            V[nx - 1, j] = V[nx - 2, j]
        elif wr == 2:
            U[nx - 1, j] = -U[nx - 2, j]
            V[nx - 1, j] = -V[nx - 2, j]
        elif wr in (3, 5):
            U[nx - 1, j] = U[nx - 2, j]
            V[nx - 1, j] = V[nx - 2, j]
        elif wr == 4:
            U[nx - 1, j] = U[1, j]
            V[nx - 1, j] = V[1, j]

    for i in range(nx):
        if wb == 1:
            U[i, 0] = U[i, 1]
            V[i, 0] = 0.0
        elif wb == 2:
            U[i, 0] = -U[i, 1]
            V[i, 0] = -V[i, 1]
        elif wb in (3, 5):
            U[i, 0] = U[i, 1]
            V[i, 0] = V[i, 1]
        elif wb == 4:
            U[i, 0] = U[i, ny - 2]
            V[i, 0] = V[i, ny - 2]

        if wt == 1:
            U[i, ny - 1] = U[i, ny - 2]
            V[i, ny - 1] = 0.0
        elif wt == 2:
            U[i, ny - 1] = -U[i, ny - 2]
            V[i, ny - 1] = -V[i, ny - 2]
        elif wt in (3, 5):
            U[i, ny - 1] = U[i, ny - 2]
            V[i, ny - 1] = V[i, ny - 2]
        elif wt == 4:
            U[i, ny - 1] = U[i, 1]
            V[i, ny - 1] = V[i, 1]

    print(f"RTINIT: NMAX = {nmax} modes included")
    print(f"RTINIT: KC = {kc:.6f}  KMAX = {kmax:.6f}")


def _dalziel_perturb(mesh: Mesh, fields: Fields, cfg: RunConfig) -> None:
    """
    Port of Fortran SUBROUTINE DALZIEL_PERTURB (lines 2470-2551).
    Dalziel et al. (1999) interface perturbation.
    Modifies F field directly. Wavelengths 4*dx to 8*dx.
    F=1 (heavy fluid) ABOVE interface, F=0 below.
    """
    F = fields.F
    xi = mesh.xi
    yj = mesh.yj
    delx = mesh.delx
    dely = mesh.dely
    im1 = mesh.im1
    jm1 = mesh.jm1

    hb = cfg.h_dist
    l = cfg.r_dist
    pi = np.pi
    twopi = 2.0 * pi

    dx = delx[1]  # Fortran DELX(2)
    sigma_p = 0.08 * dx

    # Mode range: wavelengths from 4*dx to 8*dx
    nmin = int(l / (8.0 * dx)) + 1
    nmax = int(l / (4.0 * dx))
    nmodes = nmax - nmin + 1

    if nmodes < 1:
        print("DALZIEL_PERTURB: ERROR - no modes in range")
        return
    if nmodes > 200:
        print("DALZIEL_PERTURB: ERROR - too many modes")
        nmodes = 200

    # Uniform amplitude, random phases (portable LCG)
    iseed = 123
    amp = np.empty(nmodes, dtype=float)
    phs = np.empty(nmodes, dtype=float)

    for n in range(nmodes):
        iseed = (iseed * 1103515245 + 12345) % 2147483647
        phs[n] = twopi * float(iseed) / 2147483647.0
        amp[n] = sigma_p * np.sqrt(2.0 / float(nmodes))

    print(f"DALZIEL_PERTURB: Dalziel 1999 interface perturbation")
    print(f"  dx = {dx}  sigma = {sigma_p}")
    print(f"  Modes: {nmin} to {nmax} ({nmodes} modes)")
    print(f"  Interface height HB = {hb}")

    # Modify F field
    # Fortran: I=2..IM1 → Python i=1..im1-1
    for i in range(1, im1):
        eta = 0.0
        for n in range(nmodes):
            kx = twopi * float(nmin + n) / l
            eta += amp[n] * np.cos(kx * xi[i] + phs[n])

        for j in range(1, jm1):
            ytop = yj[j] + 0.5 * dely[j]
            ybot = yj[j] - 0.5 * dely[j]

            if ytop <= hb + eta:
                # Cell entirely below interface: light fluid
                F[i, j] = 0.0
            elif ybot >= hb + eta:
                # Cell entirely above interface: heavy fluid
                F[i, j] = 1.0
            else:
                # Cell straddles interface: fractional fill
                frac = (ytop - (hb + eta)) / dely[j]
                F[i, j] = max(0.0, min(1.0, frac))

