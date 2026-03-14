# io.py
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple
import numpy as np

from .mesh import Mesh
from .fields import Fields
from .state import RunState

from .mesh import MeshInput
from .config import RunConfig

def _read_bool(token: str) -> bool:
    token = token.strip()
    if token.lower() in ("true", "t", ".true."):
        return True
    if token.lower() in ("false", "f", ".false."):
        return False
    raise ValueError(f"Cannot parse boolean from '{token}'")


def read_keyword_deck(path: str | Path) -> Tuple[MeshInput, RunConfig]:
    """
    Parse a keyword-style deck like int200.in into MeshInput + RunConfig.[file:57]
    """
    path = Path(path)
    lines = [ln.split("#", 1)[0].strip() for ln in path.read_text().splitlines()]
    # drop empty lines
    tokens = [ln for ln in lines if ln]

    def pop_line() -> str:
        if not tokens:
            raise ValueError("Unexpected end of file")
        return tokens.pop(0)

    # --- mesh block ---
    assert pop_line().lower() == "mesh"
    icyl_raw = pop_line().split()[0]
    # ICYL can be bool (true/false) or int (0/1) in input decks
    try:
        icyl = int(icyl_raw)
    except ValueError:
        icyl = 1 if _read_bool(icyl_raw) else 0
    nkx, nky = map(int, pop_line().split())

    xc = np.fromstring(pop_line(), sep=" ", dtype=float)
    assert xc.size == nkx

    nxl = np.fromstring(pop_line(), sep=" ", dtype=int)
    nxr = np.fromstring(pop_line(), sep=" ", dtype=int)
    dxmn = np.fromstring(pop_line(), sep=" ", dtype=float)
    yc = np.fromstring(pop_line(), sep=" ", dtype=float)
    nyl = np.fromstring(pop_line(), sep=" ", dtype=int)
    nyr = np.fromstring(pop_line(), sep=" ", dtype=int)
    dymn = np.fromstring(pop_line(), sep=" ", dtype=float)

    xl = np.fromstring(pop_line(), sep=" ", dtype=float)
    yl = np.fromstring(pop_line(), sep=" ", dtype=float)

    ifl, ifr, jfb, jft = map(int, pop_line().split())

    # --- properties block ---
    assert pop_line().lower() == "properties"
    nmat = int(pop_line().split()[0])
    vnu, vnuc = map(float, pop_line().split())
    rhof, rhofc = map(float, pop_line().split())
    isurf10_str, sigma_str, cangle_str = pop_line().split()
    isurf10 = _read_bool(isurf10_str)
    sigma = float(sigma_str)
    cangle = float(cangle_str)

    # --- output block ---
    assert pop_line().lower() == "output"
    ip_str, iare_str, imovy_str = pop_line().split()
    iplots = _read_bool(ip_str)
    iare = _read_bool(iare_str)
    imovy = _read_bool(imovy_str)
    prtdt = float(pop_line())
    pltst1, pltst2 = map(float, pop_line().split())
    pltdt1, pltdt2 = map(float, pop_line().split())

    # Optional prefix line (default "svof-")
    # Peek: if the next token is the "bcs" section header, use default
    next_line = pop_line()
    if next_line.lower() == "bcs":
        prefix = "svof-"
    else:
        prefix = next_line.strip()
        # Now consume the actual "bcs" header
        assert pop_line().lower() == "bcs"

    # --- bcs block ---
    iwl, iwr, iwt, iwb = map(int, pop_line().split())

    # --- particles block ---
    assert pop_line().lower() == "particles"
    npx, npy = map(int, pop_line().split())
    xpl, xpr = map(float, pop_line().split())
    ypb, ypt = map(float, pop_line().split())

    # --- obstacles block ---
    assert pop_line().lower() == "obstacles"
    nobs = int(pop_line().split()[0])
    iomin, iomax = map(int, pop_line().split())
    jomin, jomax = map(int, pop_line().split())

    # --- disturb block ---
    assert pop_line().lower() == "disturb"
    ndis = int(pop_line().split()[0])
    idis = int(pop_line().split()[0])
    vorig = float(pop_line().split()[0])
    r_dist, h_dist, ilow, ihigh = pop_line().split()
    r_dist = float(r_dist)
    h_dist = float(h_dist)
    ilow = int(ilow)
    ihigh = int(ihigh)

    # --- svof block ---
    assert pop_line().lower() == "svof"
    epsi = float(pop_line().split()[0])
    omg = float(pop_line().split()[0])
    alpha = float(pop_line().split()[0])
    gx, gy = map(float, pop_line().split())
    ui, vi = map(float, pop_line().split())
    csq = float(pop_line().split()[0])
    fmass, flength, ftime = map(float, pop_line().split())
    ishvel_str, shvel_str = pop_line().split()
    ishvel = _read_bool(ishvel_str)
    shvel = float(shvel_str)
    delt, autot = map(float, pop_line().split())
    twfin = float(pop_line().split()[0])

    mesh_input = MeshInput(
        nkx=nkx, nky=nky,
        xl=xl, xc=xc, nxl=nxl, nxr=nxr, dxmn=dxmn,
        yl=yl, yc=yc, nyl=nyl, nyr=nyr, dymn=dymn,
        iwl=iwl, iwr=iwr, iwt=iwt, iwb=iwb
    )

    run_config = RunConfig(
        icyl=icyl,
        nmat=nmat, vnu=vnu, vnuc=vnuc, rhof=rhof, rhofc=rhofc,
        isurf10=isurf10, sigma=sigma, cangle_deg=cangle,
        ifl=ifl, ifr=ifr, jfb=jfb, jft=jft,
        prefix=prefix,
        iplots=iplots, iare=iare, imovy=imovy, prtdt=prtdt,
        pltst1=pltst1, pltst2=pltst2, pltdt1=pltdt1, pltdt2=pltdt2,
        iwl=iwl, iwr=iwr, iwt=iwt, iwb=iwb,
        npx=npx, npy=npy, xpl=xpl, xpr=xpr, ypb=ypb, ypt=ypt,
        nobs=nobs, iomin=iomin, iomax=iomax, jomin=jomin, jomax=jomax,
        ndis=ndis, idis=idis, vorig=vorig,
        r_dist=r_dist, h_dist=h_dist, ilow=ilow, ihigh=ihigh,
        epsi=epsi, omg=omg, alpha=alpha,
        gx=gx, gy=gy, ui=ui, vi=vi, csq=csq,
        fmass=fmass, flength=flength, ftime=ftime,
        ishvel=ishvel, shvel=shvel,
        delt=delt, autot=autot, twfin=twfin
    )

    return mesh_input, run_config
# io.py

def _write_array_5col(fh, arr_2d, ibar):
    """Write a 2D field in Fortran order: J outer, I inner, 5 values per line, E15.8."""
    ncols = arr_2d.shape[0]   # ibar
    nrows = arr_2d.shape[1]   # jbar
    vals_per_line = 5
    for j in range(nrows):
        for k_start in range(0, ncols, vals_per_line):
            k_end = min(k_start + vals_per_line, ncols)
            line = " ".join(f"{arr_2d[i, j]:15.8E}" for i in range(k_start, k_end))
            fh.write(line + "\n")


def write_vtk_snapshot(outdir: Path, mesh: Mesh, fields: Fields,
                       state: RunState, step: int,
                       prefix: str = "svof-") -> None:
    """
    Write a VTK RECTILINEAR_GRID snapshot matching the Fortran WRITEVTK format:
      - node-based coordinates (IBAR+1 x JBAR+1) from mesh.x, mesh.y
      - CELL_DATA with IBAR*JBAR cells
      - F, P, U, V as separate SCALARS, 5 values per line in E15.8 format
      - filename: {prefix}{ibar}x{jbar}-{time_ms}.vtk
    """
    outdir = Path(outdir)
    ibar, jbar = mesh.ibar, mesh.jbar
    time_index = int(state.t * 1000.0)
    fname = outdir / f"{prefix}{ibar}x{jbar}-{time_index}.vtk"

    # Node counts (one more than cell counts)
    nx = ibar + 1
    ny = jbar + 1
    nz = 1

    # Node positions: mesh.x[1..ibar+1] are the physical cell faces
    x_nodes = mesh.x[1:1 + nx]
    y_nodes = mesh.y[1:1 + ny]

    # Interior cell data: indices 1..ibar, 1..jbar in 0-based arrays
    # (matches Fortran I=2..IBAR+1, J=2..JBAR+1)
    sl = np.s_[1:1 + ibar, 1:1 + jbar]
    F_phys = fields.F[sl]
    P_phys = fields.P[sl]
    U_phys = fields.U[sl]
    V_phys = fields.V[sl]

    ncells = ibar * jbar

    with fname.open("w") as f:
        f.write("# vtk DataFile Version 2.0\n")
        f.write("R-T Fractals\n")
        f.write("ASCII\n")
        f.write(" DATASET RECTILINEAR_GRID\n")
        f.write(f"DIMENSIONS {nx:5d} {ny:5d} {nz:5d}\n")

        # X coordinates (node positions, 5 per line)
        f.write(f"X_COORDINATES {nx:5d} float\n")
        for k_start in range(0, nx, 5):
            k_end = min(k_start + 5, nx)
            line = " ".join(f"{x_nodes[i]:15.8E}" for i in range(k_start, k_end))
            f.write(line + "\n")

        # Y coordinates
        f.write(f"Y_COORDINATES {ny:5d} float\n")
        for k_start in range(0, ny, 5):
            k_end = min(k_start + 5, ny)
            line = " ".join(f"{y_nodes[j]:15.8E}" for j in range(k_start, k_end))
            f.write(line + "\n")

        # Z coordinate
        f.write(f"Z_COORDINATES {nz:5d} float\n")
        f.write(f"{0.0:15.8E}\n")

        f.write(f"CELL_DATA {ncells:8d}\n")

        for name, data in [("F", F_phys), ("P", P_phys),
                           ("U", U_phys), ("V", V_phys)]:
            f.write(f"SCALARS {name} float  1\n")
            f.write("LOOKUP_TABLE default\n")
            _write_array_5col(f, data, ibar)

    print(f"Wrote VTK snapshot: {fname}")

