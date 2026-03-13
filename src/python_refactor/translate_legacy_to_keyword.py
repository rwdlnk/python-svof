#!/usr/bin/env python3
"""
Translate a legacy SOLA-VOF Fortran-style input (intXXX.in)
to the keyword-style deck used by the Python refactor.[file:57]
"""

import argparse
from pathlib import Path


def parse_legacy_deck(path: Path) -> dict:
    """Parse the specific legacy format shown in int200.in into a dict.[file:57]"""
    text = path.read_text().splitlines()

    # helper to find line starting with a given token
    def find_value(prefix):
        for ln in text:
            if ln.strip().startswith(prefix):
                # split at '=' and take the right-hand side
                rhs = ln.split("=", 1)[1]
                return rhs.strip()
        raise ValueError(f"Could not find line starting with '{prefix}'")

    data = {}

    # Properties, BC, materials, etc. from labeled lines
    data["IPLOT"] = int(find_value("IPLOT").split()[0])

    data["NU"] = float(find_value("NU").split()[0])
    data["NUC"] = float(find_value("NUC").split()[0])
    data["ICYL"] = int(find_value("ICYL").split()[0])
    data["EPSI"] = float(find_value("EPSI").split()[0])
    data["GX"] = float(find_value("GX").split()[0])
    data["GY"] = float(find_value("GY").split()[0])
    data["UI"] = float(find_value("UI").split()[0])
    data["VI"] = float(find_value("VI").split()[0])
    data["VELMX"] = float(find_value("VELMX").split()[0])
    data["IMOVY"] = int(find_value("IMOVY").split()[0])
    data["OMG"] = float(find_value("OMG").split()[0])
    data["ALPHA"] = float(find_value("ALPHA").split()[0])
    data["WL"] = int(find_value("WL").split()[0])
    data["WR"] = int(find_value("WR").split()[0])
    data["WT"] = int(find_value("WT").split()[0])
    data["WB"] = int(find_value("WB").split()[0])
    data["CSQ"] = float(find_value("CSQ").split()[0])
    data["AUTOT"] = float(find_value("AUTOT").split()[0])
    data["ISYMPLT"] = int(find_value("ISYMPLT").split()[0])
    data["ISURF10"] = int(find_value("ISURF10").split()[0])
    data["SIGMA"] = float(find_value("SIGMA").split()[0])
    data["CANGLE"] = float(find_value("CANGLE").split()[0])
    data["NMAT"] = int(find_value("NMAT").split()[0])
    data["RHOF"] = float(find_value("RHOF").split()[0])
    data["RHOFC"] = float(find_value("RHOFC").split()[0])
    data["FLHT"] = float(find_value("FLHT").split()[0])
    data["XPL"] = float(find_value("XPL").split()[0])
    data["YPB"] = float(find_value("YPB").split()[0])
    data["XPR"] = float(find_value("XPR").split()[0])
    data["YPT"] = float(find_value("YPT").split()[0])
    data["NPX"] = int(find_value("NPX").split()[0])
    data["NPY"] = int(find_value("NPY").split()[0])
    data["NOBS"] = int(find_value("NOBS").split()[0])
    data["NDIS"] = int(find_value("NDIS").split()[0])
    data["IDIS"] = int(find_value("IDIS").split()[0])
    data["VORIG"] = float(find_value("VORIG").split()[0])
    data["ISHVEL"] = int(find_value("ISHVEL").split()[0])
    data["SHVEL"] = float(find_value("SHVEL").split()[0])

    data["DELT"] = float(find_value("DELT").split()[0])
    data["TWFIN"] = float(find_value("TWFIN").split()[0])
    data["PRTDT"] = float(find_value("PRTDT").split()[0])
    data["PLTST1"] = float(find_value("pltst1").split()[0])
    data["PLTST2"] = float(find_value("pltst2").split()[0])
    data["PLTDT1"] = float(find_value("PLTDT1").split()[0])
    data["PLTDT2"] = float(find_value("PLTDT2").split()[0])

    # Mesh block at the bottom: NKX, NKY, NONUNIF, XC,NXL,NXR,DXMN, YC,NYL,NYR,DYMN, XL, YL, IFL,IFR,JFB,JFT, R,H,ILOW,IHIGH[f57]
    # Find the index of "-NKX!-NKY!NONU!"
    idx = None
    for k, ln in enumerate(text):
        if "-NKX!-NKY!NONU!" in ln:
            idx = k
            break
    if idx is None:
        raise ValueError("Could not find '-NKX!-NKY!NONU!' line in legacy deck")

    nkx_nky_line = text[idx + 1]
    data["NKX"], data["NKY"], data["NONUNIF"] = map(int, nkx_nky_line.split())

    xc_line = text[idx + 3]
    yc_line = text[idx + 5]
    data["XC"] = float(xc_line.split()[0])
    data["NXL"] = int(xc_line.split()[1])
    data["NXR"] = int(xc_line.split()[2])
    data["DXMN"] = float(xc_line.split()[3])

    data["YC"] = float(yc_line.split()[0])
    data["NYL"] = int(yc_line.split()[1])
    data["NYR"] = int(yc_line.split()[2])
    data["DYMN"] = float(yc_line.split()[3])

    # XL and YL lines
    # After a comment line of dashes, the next line has XL values
    xl_line = text[idx + 9]
    yl_line = text[idx + 12]
    data["XL_vals"] = xl_line.split()
    data["YL_vals"] = yl_line.split()

    # IFL, IFR, JFB, JFT
    ifl_line = text[idx + 15]
    data["IFL"], data["IFR"], data["JFB"], data["JFT"] = map(int, ifl_line.split())

    # Disturbance R,H,ILOW,IHIGH
    disturb_line = text[idx + 19]
    parts = disturb_line.split()
    data["R_dist"] = float(parts[0])
    data["H_dist"] = float(parts[1])
    data["ILOW"] = int(parts[2])
    data["IHIGH"] = int(parts[3])

    return data


def write_keyword_deck(out_path: Path, d: dict) -> None:
    """Write keyword-style deck using the parsed data."""
    with out_path.open("w") as f:
        f.write("# Converted from legacy SOLA-VOF deck\n")

        # mesh block
        f.write("mesh\n")
        f.write("False # icyl\n")
        f.write(f"{d['NKX']} {d['NKY']} # NKX NKY\n")
        f.write(f"{d['XC']}  # XC (1D array, NKX long)\n")
        f.write(f"{d['NXL']}  # NXL (1D array, NKX long)\n")
        f.write(f"{d['NXR']}  # NXR (1D array, NKX long)\n")
        f.write(f"{d['DXMN']} # DXmin (1D array, NKX long)\n")
        f.write(f"{d['YC']}  # YC (1D array, NKY long)\n")
        f.write(f"{d['NYL']}  # NYL (1D array, NKY long)\n")
        f.write(f"{d['NYR']}  # NYR (1D array, NKY long)\n")
        f.write(f"{d['DYMN']} # DYmin (1D array, NKY long)\n")
        f.write(" ".join(d['XL_vals']) + " # XL (1D array, NKX+1 long)\n")
        f.write(" ".join(d['YL_vals']) + " # YL (1D array, NKY+1 long)\n")
        f.write(f"{d['IFL']} {d['IFR']} {d['JFB']} {d['JFT']} # IFL IFR JFB JFT\n")
        f.write("#\n#\n")

        # properties
        f.write("properties\n")
        f.write(f"{d['NMAT']} # NMAT\n")
        f.write(f"{d['NU']} {d['NUC']} # NU NUC\n")
        f.write(f"{d['RHOF']} {d['RHOFC']} # RHOF RHOFC\n")
        isurf10_bool = "True" if d["ISURF10"] != 0 else "False"
        f.write(f"{isurf10_bool} {d['SIGMA']} {d['CANGLE']} # ISURF10 sigma cangle\n")
        f.write("#\n")

        # output
        f.write("output\n")
        iplots = "True" if d["IPLOT"] != 0 else "False"
        f.write(f"{iplots} False False # IPLOTS IARE IMOVY\n")
        f.write(f"{d['PRTDT']} # PRTDT\n")
        f.write(f"{d['PLTST1']} {d['PLTST2']} # PLTST1 PLTST2\n")
        f.write(f"{d['PLTDT1']} {d['PLTDT2']} # PLTDT1 PLTDT2\n")
        f.write("#\n")

        # bcs
        f.write("bcs\n")
        f.write(f"{d['WL']} {d['WR']} {d['WT']} {d['WB']} # WL WR WT WB\n")
        f.write("#\n")

        # particles
        f.write("particles\n")
        f.write(f"{d['NPX']} {d['NPY']} # NPX NPY\n")
        f.write(f"{d['XPL']} {d['XPR']} # XPL XPR\n")
        f.write(f"{d['YPB']} {d['YPT']} # YPB YPT\n")
        f.write("#\n")

        # obstacles
        f.write("obstacles\n")
        f.write(f"{d['NOBS']} # NOBS\n")
        f.write("1 1 # IOMIN IOMAX (placeholder)\n")
        f.write("1 1 # JOMIN JOMAX (placeholder)\n")
        f.write("#\n")

        # disturb
        f.write("disturb\n")
        f.write(f"{d['NDIS']} # NDIS\n")
        f.write(f"{d['IDIS']} # IDIS\n")
        f.write(f"{d['VORIG']} # VORIG\n")
        f.write(f"{d['R_dist']} {d['H_dist']} {d['ILOW']} {d['IHIGH']} # R H ILOW IHIGH\n")
        f.write("#\n")

        # svof
        f.write("svof\n")
        f.write(f"{d['EPSI']} # EPSI\n")
        f.write(f"{d['OMG']} # OMG\n")
        f.write(f"{d['ALPHA']} # ALPHA\n")
        f.write(f"{d['GX']} {d['GY']} # Gx Gy\n")
        f.write(f"{d['UI']} {d['VI']} # UI VI\n")
        f.write(f"{d['CSQ']} # CSQ\n")
        f.write("1. 1. 1. # Fmass Flength Ftime\n")
        ishvel_bool = "True" if d["ISHVEL"] != 0 else "False"
        f.write(f"{ishvel_bool} {d['SHVEL']} # ISHvel SHvel\n")
        f.write(f"{d['DELT']} {d['AUTOT']} # DELT AUTOT\n")
        f.write(f"{d['TWFIN']} # TWfin\n")

def main():
    ap = argparse.ArgumentParser(description="Translate legacy SOLA-VOF deck to keyword format")
    ap.add_argument("-i", "--input", required=True, help="Legacy intXXX.in path")
    ap.add_argument("-o", "--output", required=True, help="Output keyword deck path")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)

    data = parse_legacy_deck(in_path)
    write_keyword_deck(out_path, data)
    print(f"Wrote keyword deck to {out_path}")

if __name__ == "__main__":
    main()

