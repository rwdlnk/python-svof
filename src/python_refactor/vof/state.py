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

    @classmethod
    def from_config(cls, cfg: RunConfig) -> "RunState":
        return cls(
            icyle=0,
            iter=0,
            nocon=0,
            nflgc=0,
            nmat=cfg.nmat,
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
        )

