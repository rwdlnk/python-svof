# state.py
from dataclasses import dataclass
from .config import RunConfig


@dataclass
class RunState:
    # integer counters / sizes (subset; we can extend as needed)
    icyle: int
    iter: int
    nocon: int
    nflgc: int
    nmat: int

    # time-like quantities
    t: float
    delt: float
    twprt: float
    twplt: float
    twfin: float

    # control flags
    autot: float
    flg: float
    flgc: float
    fnoc: float

    # geometry
    cyl: float          # 0.0=planar, 1.0=cylindrical

    # physical parameters (copied from RunConfig for convenience)
    gx: float
    gy: float
    vnu: float
    vnuc: float
    rhof: float
    rhofc: float
    sigma: float
    cangle_deg: float
    csq: float
    epsi: float
    omg: float
    alpha: float
    isurf10: bool
    ishvel: bool
    shvel: float
    vchgt: float

    # Stability limits (set by compute_beta)
    dtvis: float     # viscous stability limit
    dtsft: float     # surface tension stability limit
    rdtexp: float    # acoustic CFL rate

    @classmethod
    def from_config(cls, cfg: RunConfig) -> "RunState":
        return cls(
            icyle=0,
            iter=0,
            nocon=0,
            nflgc=0,
            nmat=cfg.nmat,
            cyl=float(cfg.icyl),
            t=0.0,
            delt=cfg.delt,
            twprt=cfg.prtdt,
            twplt=cfg.pltst1,
            twfin=cfg.twfin,
            autot=cfg.autot,
            flg=1.0,
            flgc=0.0,
            fnoc=0.0,
            gx=cfg.gx,
            gy=cfg.gy,
            vnu=cfg.vnu,
            vnuc=cfg.vnuc,
            rhof=cfg.rhof,
            rhofc=cfg.rhofc,
            sigma=cfg.sigma,
            cangle_deg=cfg.cangle_deg,
            csq=cfg.csq,
            epsi=cfg.epsi,
            omg=cfg.omg,
            alpha=cfg.alpha,
            isurf10=cfg.isurf10,
            ishvel=cfg.ishvel,
            shvel=cfg.shvel,
            vchgt=0.0,
            dtvis=1.0e10,
            dtsft=1.0e10,
            rdtexp=1.0e10,
        )

