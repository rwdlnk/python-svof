# config.py
from dataclasses import dataclass
from typing import Optional


@dataclass
class RunConfig:
    # geometry
    icyl: int           # 0=planar, 1=cylindrical

    # materials / properties
    nmat: int
    vnu: float
    vnuc: float
    rhof: float
    rhofc: float
    isurf10: bool
    sigma: float
    cangle_deg: float

    # output / timing
    prefix: str         # output filename prefix (default "svof-")
    iplots: bool
    iare: bool
    imovy: bool
    prtdt: float
    pltst1: float
    pltst2: float
    pltdt1: float
    pltdt2: float

    # fluid region (F=1) index bounds (1-based Fortran indices)
    ifl: int
    ifr: int
    jfb: int
    jft: int

    # boundary conditions
    iwl: int
    iwr: int
    iwt: int
    iwb: int

    # particles
    npx: int
    npy: int
    xpl: float
    xpr: float
    ypb: float
    ypt: float

    # obstacles (for now, one block; can generalize)
    nobs: int
    iomin: int
    iomax: int
    jomin: int
    jomax: int

    # disturbances
    ndis: int
    idis: int
    vorig: float
    r_dist: float
    h_dist: float
    ilow: int
    ihigh: int

    # SVOF / solver controls
    epsi: float
    omg: float
    alpha: float
    gx: float
    gy: float
    ui: float
    vi: float
    csq: float
    fmass: float
    flength: float
    ftime: float
    ishvel: bool
    shvel: float
    delt: float
    autot: float
    twfin: float

