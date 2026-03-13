# fields.py
from dataclasses import dataclass
import numpy as np
from .mesh import Mesh


@dataclass
class Fields:
    # time-n fields
    UN: np.ndarray
    VN: np.ndarray
    PN: np.ndarray
    FN: np.ndarray

    # current fields
    U: np.ndarray
    V: np.ndarray
    P: np.ndarray
    F: np.ndarray
    PETA: np.ndarray
    BETA: np.ndarray
    DTANTH: np.ndarray
    PS: np.ndarray
    NF: np.ndarray     # integer flags

    Fbar: np.ndarray   # e.g., shape (jbar-1,) or similar

    @classmethod
    def allocate(cls, mesh: Mesh) -> "Fields":
        """
        Allocate field arrays sized by the mesh (including ghost cells).[file:41]
        """
        imax, jmax = mesh.imax, mesh.jmax
        shape = (imax, jmax)

        UN = np.zeros(shape, dtype=float)
        VN = np.zeros(shape, dtype=float)
        PN = np.zeros(shape, dtype=float)
        FN = np.zeros(shape, dtype=float)

        U = np.zeros(shape, dtype=float)
        V = np.zeros(shape, dtype=float)
        P = np.zeros(shape, dtype=float)
        F = np.zeros(shape, dtype=float)
        PETA = np.ones(shape, dtype=float)
        BETA = np.zeros(shape, dtype=float)
        DTANTH = np.zeros(shape, dtype=float)
        PS = np.zeros(shape, dtype=float)
        NF = np.zeros(shape, dtype=int)

        # Fbar corresponds to 1D array over J (like NJF in Fortran)
        # Here we take it as jbar-1; you can tweak once we wire PRT, etc.
        njf = mesh.jbar - 1 if mesh.jbar > 1 else 1
        Fbar = np.zeros(njf, dtype=float)

        return cls(
            UN=UN, VN=VN, PN=PN, FN=FN,
            U=U, V=V, P=P, F=F,
            PETA=PETA, BETA=BETA, DTANTH=DTANTH, PS=PS, NF=NF,
            Fbar=Fbar
        )

