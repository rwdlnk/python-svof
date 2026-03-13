# mesh.py
from dataclasses import dataclass
from typing import Optional
import numpy as np

@dataclass
class Mesh:
    # 1D coordinate arrays including ghost cells
    x: np.ndarray       # cell-face x
    xi: np.ndarray      # cell-center x
    delx: np.ndarray
    rdx: np.ndarray
    rx: np.ndarray

    y: np.ndarray       # cell-face y
    yj: np.ndarray      # cell-center y
    dely: np.ndarray
    rdy: np.ndarray
    ryj: np.ndarray

    # logical mesh parameters
    ibar: int           # number of physical cells in x
    jbar: int           # number of physical cells in y
    imax: int           # total cells in x including ghosts
    jmax: int           # total cells in y including ghosts
    im1: int
    jm1: int
    im2: int
    jm2: int

    # bookkeeping for nonuniform blocks (Fortran NKX, NKY, etc.)
    nkx: int
    nky: int
    xl: np.ndarray
    xc: np.ndarray
    dxmn: np.ndarray
    nxl: np.ndarray
    nxr: np.ndarray
    yl: np.ndarray
    yc: np.ndarray
    dymn: np.ndarray
    nyl: np.ndarray
    nyr: np.ndarray

    # periodicity flags (derived from IWR, IWT)
    periodic_x: bool
    periodic_y: bool


@dataclass
class MeshInput:
    """
    Minimal container for the mesh-related inputs that in Fortran are read
    before MESHSET is called: NKX, NKY, NONUNIF, XL, XC, NXL, NXR, DXMN, etc.
    """
    nkx: int
    nky: int
    xl: np.ndarray      # length nkx+1
    xc: np.ndarray      # length nkx
    nxl: np.ndarray     # length nkx
    nxr: np.ndarray     # length nkx
    dxmn: np.ndarray    # length nkx
    yl: np.ndarray      # length nky+1
    yc: np.ndarray      # length nky
    nyl: np.ndarray     # length nky
    nyr: np.ndarray     # length nky
    dymn: np.ndarray    # length nky
    iwl: int            # left wall type (1–5)
    iwr: int            # right wall type (1–5)
    iwt: int            # top wall type (1–5)
    iwb: int            # bottom wall type (1–5)


def meshset(mesh_in: MeshInput) -> Mesh:
    """
    Python version of Fortran MESHSET, returning a Mesh object with ghost zones.

    It reproduces the piecewise-nonuniform construction in x and y,
    then adds 2 ghost cells on each nonperiodic side and an extra
    ghost on periodic sides (matching IWR/IWT=4 behavior).[file:40]
    """
    nkx = mesh_in.nkx
    nky = mesh_in.nky

    # --- Build physical-face coordinates X, Y exactly as in MESHSET ---

    # First build physical 1D X faces (no ghosts yet)
    # Count number of physical faces in x
    # Fortran logic walks through each block K, adding NXL(K)+NXR(K) faces
    numx = 1
    for k in range(nkx):
        numx += mesh_in.nxl[k] + mesh_in.nxr[k]
    x_phys = np.empty(numx, dtype=float)

    i = 0
    x_phys[i] = mesh_in.xl[0]
    for k in range(nkx):
        dxml = (mesh_in.xc[k] - mesh_in.xl[k]) / mesh_in.nxl[k]
        dxmr = (mesh_in.xl[k + 1] - mesh_in.xc[k]) / mesh_in.nxr[k]
        dxmn1 = mesh_in.dxmn[k]

        # left side of block
        nt = mesh_in.nxl[k]
        tn = float(nt)
        tn = max(tn, 1.0 + 1.0e-6)  # EM6
        dxmn = min(dxmn1, dxml)
        cmc = (mesh_in.xc[k] - mesh_in.xl[k] - tn * dxmn) * tn / (tn - 1.0)
        if nt == 1:
            cmc = 0.0
        bmc = mesh_in.xc[k] - mesh_in.xl[k] - cmc
        for l in range(1, nt + 1):
            i += 1
            rln = (float(l) - tn) / tn
            x_phys[i] = mesh_in.xc[k] + bmc * rln - cmc * rln * rln

        # right side of block
        nt = mesh_in.nxr[k]
        tn = float(nt)
        tn = max(tn, 1.0 + 1.0e-6)
        dxmn = min(dxmn1, dxmr)
        cmc = (mesh_in.xl[k + 1] - mesh_in.xc[k] - tn * dxmn) * tn / (tn - 1.0)
        if nt == 1:
            cmc = 0.0
        bmc = mesh_in.xl[k + 1] - mesh_in.xc[k] - cmc
        for l in range(1, nt + 1):
            i += 1
            rln = float(l) / tn
            x_phys[i] = mesh_in.xc[k] + bmc * rln + cmc * rln * rln

    # Fortran may add an extra cell for IWR=4 (periodic); we handle via ghosts
    numx_phys = numx  # number of physical faces

    # Build physical 1D Y faces similarly
    numy = 1
    for k in range(nky):
        numy += mesh_in.nyl[k] + mesh_in.nyr[k]
    y_phys = np.empty(numy, dtype=float)

    j = 0
    y_phys[j] = mesh_in.yl[0]
    for k in range(nky):
        dyml = (mesh_in.yc[k] - mesh_in.yl[k]) / mesh_in.nyl[k]
        dymr = (mesh_in.yl[k + 1] - mesh_in.yc[k]) / mesh_in.nyr[k]
        dymn1 = mesh_in.dymn[k]

        # bottom half of block
        nt = mesh_in.nyl[k]
        tn = float(nt)
        tn = max(tn, 1.0 + 1.0e-6)
        dymn = min(dymn1, dyml)
        cmc = (mesh_in.yc[k] - mesh_in.yl[k] - tn * dymn) * tn / (tn - 1.0)
        if nt == 1:
            cmc = 0.0
        bmc = mesh_in.yc[k] - mesh_in.yl[k] - cmc
        for l in range(1, nt + 1):
            j += 1
            rln = (float(l) - tn) / tn
            y_phys[j] = mesh_in.yc[k] + bmc * rln - cmc * rln * rln

        # top half of block
        nt = mesh_in.nyr[k]
        tn = float(nt)
        tn = max(tn, 1.0 + 1.0e-6)
        dymn = min(dymn1, dymr)
        cmc = (mesh_in.yl[k + 1] - mesh_in.yc[k] - tn * dymn) * tn / (tn - 1.0)
        if nt == 1:
            cmc = 0.0
        bmc = mesh_in.yl[k + 1] - mesh_in.yc[k] - cmc
        for l in range(1, nt + 1):
            j += 1
            rln = float(l) / tn
            y_phys[j] = mesh_in.yc[k] + bmc * rln + cmc * rln * rln

    numy_phys = numy

    # ibar/jbar: number of internal cells (like IBAR, JBAR)
    ibar = numx_phys - 1
    jbar = numy_phys - 1

    # Periodicity from wall types (IWR/IWT=4 => periodic)
    periodic_x = (mesh_in.iwr == 4)
    periodic_y = (mesh_in.iwt == 4)

    # Ghost-cell counts
    nghost_x = 2 + (1 if periodic_x else 0)
    nghost_y = 2 + (1 if periodic_y else 0)

    # Total cells including ghosts
    imax = ibar + nghost_x
    jmax = jbar + nghost_y

    # Build x with ghosts: prepend/append from physical faces
    # For simplicity, make ghosts mirror spacing at boundaries;
    # periodic will be enforced later in BC.
    dx_left = x_phys[1] - x_phys[0]
    dx_right = x_phys[-1] - x_phys[-2]
    x = np.empty(imax, dtype=float)
    # interior (physical) region
    x[1:1 + numx_phys] = x_phys
    # left ghost(s)
    x[0] = x[1] - dx_left
    if nghost_x == 3:
        # one more ghost on the right for periodic case
        x[1 + numx_phys] = x[1 + numx_phys - 1] + dx_right
    # rightmost ghost (if nonperiodic, this is the only right ghost)
    x[-1] = x[-2] + dx_right

    # Same for y
    dy_bottom = y_phys[1] - y_phys[0]
    dy_top = y_phys[-1] - y_phys[-2]
    y = np.empty(jmax, dtype=float)
    y[1:1 + numy_phys] = y_phys
    y[0] = y[1] - dy_bottom
    if nghost_y == 3:
        y[1 + numy_phys] = y[1 + numy_phys - 1] + dy_top
    y[-1] = y[-2] + dy_top

    # Derived arrays: centers, spacings, reciprocals
    xi = 0.5 * (x[:-1] + x[1:])
    yj = 0.5 * (y[:-1] + y[1:])
    delx = np.diff(x)
    dely = np.diff(y)
    rdx = 1.0 / delx
    rdy = 1.0 / dely

    rx = np.zeros_like(x)
    nonzero_x = np.abs(x) > 0.0
    rx[nonzero_x] = 1.0 / x[nonzero_x]

    ryj = np.zeros_like(yj)
    nonzero_yj = np.abs(yj) > 0.0
    ryj[nonzero_yj] = 1.0 / yj[nonzero_yj]

    im1 = imax - 1
    jm1 = jmax - 1
    im2 = imax - 2
    jm2 = jmax - 2

    return Mesh(
        x=x, xi=xi, delx=delx, rdx=rdx, rx=rx,
        y=y, yj=yj, dely=dely, rdy=rdy, ryj=ryj,
        ibar=ibar, jbar=jbar, imax=imax, jmax=jmax,
        im1=im1, jm1=jm1, im2=im2, jm2=jm2,
        nkx=nkx, nky=nky,
        xl=mesh_in.xl, xc=mesh_in.xc, dxmn=mesh_in.dxmn,
        nxl=mesh_in.nxl, nxr=mesh_in.nxr,
        yl=mesh_in.yl, yc=mesh_in.yc, dymn=mesh_in.dymn,
        nyl=mesh_in.nyl, nyr=mesh_in.nyr,
        periodic_x=periodic_x, periodic_y=periodic_y
    )

