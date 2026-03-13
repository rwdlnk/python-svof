#!/usr/bin/env python3
"""
Read SOLA-VOF or VOF (Python refactor) VTK rectilinear grid files and
compute velocity maxima.

Two VTK formats are supported:
  sola  — CELL_DATA with separate SCALARS U, V (DIMENSIONS = nodes)
  vof   — POINT_DATA with VECTORS velocity (DIMENSIONS = cells)

Usage:
  python3 vmax.py sola file.vtk
  python3 vmax.py sola RT160x200-*.vtk
  python3 vmax.py vof  vof_00200.vtk
  python3 vmax.py vof  /path/to/case/vof_*.vtk
"""

import sys
import re
import numpy as np


def read_sola_vtk(vtk_file):
    """Read SOLA-VOF VTK rectilinear grid file.

    Format: CELL_DATA, DIMENSIONS = node counts (nx+1, ny+1, 1),
    separate SCALARS for F, P, U, V.

    Returns:
        fields: dict of field_name -> 2D array (ny, nx)
        nx, ny: cell counts
    """
    with open(vtk_file, 'r') as f:
        content = f.read()

    m = re.search(r'DIMENSIONS\s+(\d+)\s+(\d+)\s+(\d+)', content)
    if not m:
        raise ValueError(f"Could not parse DIMENSIONS in {vtk_file}")
    nx = int(m.group(1)) - 1  # node count - 1 = cell count
    ny = int(m.group(2)) - 1

    fields = {}
    for m in re.finditer(
            r'SCALARS\s+(\w+)\s+\w+\s+\d+\s*\n'
            r'LOOKUP_TABLE\s+\w+\s*\n'
            r'([\s\S]*?)(?=SCALARS|\Z)', content):
        name = m.group(1)
        vals = np.array([float(v) for v in m.group(2).split()])
        if len(vals) == nx * ny:
            fields[name] = vals.reshape(ny, nx)

    return fields, nx, ny


def read_vof_vtk(vtk_file):
    """Read VOF (Python refactor) VTK rectilinear grid file.

    Format: POINT_DATA, DIMENSIONS = cell counts (nx, ny, 1),
    VECTORS velocity float (3-component: u, v, w).

    Returns:
        fields: dict with 'U' and 'V' -> 2D array (ny, nx)
        nx, ny: cell counts
    """
    with open(vtk_file, 'r') as f:
        content = f.read()

    m = re.search(r'DIMENSIONS\s+(\d+)\s+(\d+)\s+(\d+)', content)
    if not m:
        raise ValueError(f"Could not parse DIMENSIONS in {vtk_file}")
    nx = int(m.group(1))  # dimensions ARE cell counts
    ny = int(m.group(2))

    fields = {}

    # Parse VECTORS velocity
    vm = re.search(r'VECTORS\s+velocity\s+\w+\s*\n([\s\S]*?)(?=SCALARS|VECTORS|\Z)',
                   content)
    if vm:
        vals = np.array([float(v) for v in vm.group(1).split()])
        n_pts = nx * ny
        if len(vals) == 3 * n_pts:
            vecs = vals.reshape(n_pts, 3)
            fields['U'] = vecs[:, 0].reshape(ny, nx)
            fields['V'] = vecs[:, 1].reshape(ny, nx)
            fields['W'] = vecs[:, 2].reshape(ny, nx)

    # Also parse any SCALARS (e.g. F, P)
    for sm in re.finditer(
            r'SCALARS\s+(\w+)\s+\w+\s+\d+\s*\n'
            r'LOOKUP_TABLE\s+\w+\s*\n'
            r'([\s\S]*?)(?=SCALARS|VECTORS|\Z)', content):
        name = sm.group(1)
        vals = np.array([float(v) for v in sm.group(2).split()])
        if len(vals) == nx * ny:
            fields[name] = vals.reshape(ny, nx)

    return fields, nx, ny


def compute_vmax(fields):
    """Compute velocity maxima from U, V fields.

    Returns:
        (max_abs_u, max_abs_v, max_speed)
    """
    U = fields['U']
    V = fields['V']
    speed = np.sqrt(U**2 + V**2)
    return float(np.max(np.abs(U))), float(np.max(np.abs(V))), float(np.max(speed))


def main():
    if len(sys.argv) < 3 or sys.argv[1] not in ('sola', 'vof'):
        print("Usage: python3 vmax.py {sola|vof} file.vtk [file2.vtk ...]")
        print("  sola  — SOLA-VOF format (CELL_DATA, separate U/V scalars)")
        print("  vof   — VOF Python refactor (POINT_DATA, VECTORS velocity)")
        sys.exit(1)

    fmt = sys.argv[1]
    reader = read_sola_vtk if fmt == 'sola' else read_vof_vtk
    files = sys.argv[2:]
    multi = len(files) > 1

    for vtk_file in files:
        try:
            fields, nx, ny = reader(vtk_file)
        except Exception as e:
            print(f"ERROR reading {vtk_file}: {e}")
            continue

        if 'U' not in fields or 'V' not in fields:
            print(f"ERROR: {vtk_file} missing U or V field "
                  f"(found: {list(fields.keys())})")
            continue

        max_u, max_v, max_speed = compute_vmax(fields)

        if multi:
            # Extract timestep index from filename
            m = re.search(r'[\-_](\d+)\.vtk$', vtk_file)
            idx = m.group(1) if m else '?'
            print(f"[{idx:>5s}]  max|U|={max_u:.6e}  "
                  f"max|V|={max_v:.6e}  max|speed|={max_speed:.6e}")
        else:
            print(f"File: {vtk_file}")
            print(f"Format: {fmt}")
            print(f"Grid: {nx}x{ny} = {nx*ny} cells")
            print(f"max |U|     = {max_u:.6e} m/s")
            print(f"max |V|     = {max_v:.6e} m/s")
            print(f"max |speed| = {max_speed:.6e} m/s")
            print(f"U range: [{fields['U'].min():.6e}, {fields['U'].max():.6e}]")
            print(f"V range: [{fields['V'].min():.6e}, {fields['V'].max():.6e}]")


if __name__ == '__main__':
    main()
