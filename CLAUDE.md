# python-svof — Claude Development Notes

## Project Overview

Python refactor of the SOLA-VOF (Volume-Of-Fluid) Fortran CFD code (`slimMaster.f`).
Simulates two-phase incompressible flow with interface tracking on structured 2D meshes.
Licensed GPLv3.

## Directory Structure

```
src/slimMaster/          # Original Fortran code
  slimMaster.f           # Full SOLA-VOF implementation (~3200 lines)
  dint.h.160x200         # Header sized for 160x200 mesh (symlink as dint.h)
  compit                 # Compile script: gfortran -std=legacy -fdec -O4

src/python_refactor/     # Python port
  main.py                # CLI driver: time loop, plot timing, VTK output
  vof/
    config.py            # RunConfig dataclass (all input parameters)
    state.py             # RunState dataclass (runtime counters, physics params)
    mesh.py              # Structured non-uniform 2D mesh generation (meshset)
    fields.py            # Field arrays: U,V,P,F,UN,VN,PN,FN,BETA,PETA,etc.
    initcond.py          # Initial conditions: F setup, hydrostatic P, disturbances
    solver.py            # tilde_step (vectorized), pressure_iteration (SOLA)
    bc.py                # Boundary conditions (5 wall types)
    io.py                # Keyword input parser + VTK writer (Fortran-compatible)
  vmax.py                # VTK analysis tool (reads both Fortran and Python output)
  translate_legacy_to_keyword.py  # Converts old Fortran input format to keyword format

test/
  svof/                  # Fortran reference input/output
  new/                   # Python keyword-format input
  compare_fortran/       # Fortran comparison runs
  compare_python/        # Python comparison runs
```

## Key Decisions and Changes Made

### VTK Output (io.py)
Rewrote `write_vtk_snapshot` to match Fortran WRITEVTK exactly:
- CELL_DATA (not POINT_DATA) with IBAR*JBAR cells
- Node-based coordinates (IBAR+1 x JBAR+1) from mesh.x, mesh.y
- F, P, U, V as four separate SCALARS, 5 values per line in E15.8 format
- Filename: `RT{ibar}x{jbar}-{time_ms}.vtk` (same as Fortran)
- This ensures all existing analysis tools (vmax.py, ParaView) work on both outputs

### Plot Timing (main.py)
Replaced step-based `--vtk-every` with Fortran-style time-based plot control:
- Uses pltst1/pltdt1/pltst2/pltdt2 from input deck
- Writes VTK at t=0, then at pltdt1 intervals, switching to pltdt2 after pltst2

### Time-N Array Copy (main.py)
Added the critical UN=U, VN=V, PN=P, FN=F copy at the start of each cycle.
The Fortran does this at line ~204. Without it, tilde_step reads stale zeros from
UN/VN and produces wrong results.

### Hydrostatic Pressure Init (initcond.py)
Added hydrostatic pressure initialization matching Fortran SETUP (line ~2698):
- Integrates P from top down: P[i,jj] = P[i,jj+1] - gy * rhoya
- Fills ghost cell pressure
- Without this, gravity is unbalanced and the solver blows up immediately

### SOLA Pressure Iteration (solver.py)
Completely replaced the crude Poisson solver with a faithful port of Fortran PRESSIT:
- Computes divergence DIJ at each cell
- Corrects P by DELP = -DIJ (for incompressible, PETA=1, RCSQ=0)
- Immediately updates neighboring U,V using DPTC = 2*DELT*DELP
- Iterates until max|DIJ| < EPSI or 1000 iterations
- Calls BC after each sweep
- This is a Gauss-Seidel-style simultaneous P/V correction (the SOLA method)

### Vectorized tilde_step (solver.py)
Rewrote the cell-by-cell Python loops as NumPy array operations.
Gives ~100x speedup over the original scalar loops.
Covers advection (WUDS/alpha scheme), viscosity, pressure gradient, and gravity.

### PLANAR Disturbance Fix (initcond.py + slimMaster.f)
Fixed three bugs in SUBROUTINE PLANAR:
1. `V(II,JJ) = VSUM` → `V(II,JJ) = V(II,JJ) + VSUM` (accumulate modes, not overwrite)
2. Random coefficients changed from [0,1) to [-1,+1] per Youngs (1984) formulation
3. VORIG now scales the total sum (not each term), so it controls max amplitude

## Known Issues / Incomplete

- **pressure_iteration is scalar loops** — the SOLA Gauss-Seidel sweep is sequential;
  needs Numba JIT or red-black vectorization for performance
- **bc.py indexing** — some ghost cell copies use `P[imax-1,j] = P[im1,j]` where
  imax-1 == im1 in 0-based indexing (copies to self). Works because hydrostatic init
  fills ghosts separately, but should be audited
- **No VOF advection** — F field is never updated during time-stepping
- **No surface tension** — sigma parsed but unused
- **No obstacles** — BETA mask allocated but never populated
- **dely sizing** — Python dely has imax-1 entries (np.diff), Fortran DELY has IBAR2;
  top ghost cell dely is clamped with min() in hydrostatic init

## Build and Run

### Fortran
```bash
cd src/slimMaster
ln -sf dint.h.160x200 dint.h
bash compit
cd ../../test/compare_fortran
../../src/slimMaster/slimMaster
```

### Python
```bash
cd src/python_refactor
python3 main.py -i ../../test/new/int200.in -o ../../test/compare_python
```

### Compare
```bash
cd src/python_refactor
python3 vmax.py sola ../../test/compare_fortran/RT160x200-*.vtk
python3 vmax.py sola ../../test/compare_python/RT160x200-*.vtk
```

## Dependencies
- Python 3.10+
- NumPy
- Numba (for JIT-accelerated pressure solver)
- gfortran (for Fortran reference code)
