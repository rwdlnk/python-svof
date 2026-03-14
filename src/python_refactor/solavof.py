# main.py
import argparse
from pathlib import Path

from vof.mesh import meshset
from vof.fields import Fields
from vof.io import read_keyword_deck, write_vtk_snapshot
from vof.state import RunState
from vof.bc import apply_boundary_conditions, free_surface_bc
from vof.initcond import setup_initial_state
from vof.solver import (tilde_step, pressure_iteration, vfconv, petacal,
                        compute_beta, deltadj, tms10)


def parse_args():
    parser = argparse.ArgumentParser(description="SOLA-VOF Python driver")
    parser.add_argument(
        "-i", "--input",
        type=str,
        required=True,
        help="Path to keyword-style input deck (e.g. int200.in)"
    )
    parser.add_argument(
        "-o", "--outdir",
        type=str,
        default=".",
        help="Output directory (default: current directory)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    deck_path = Path(args.input).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Read deck and build mesh + config
    mesh_input, cfg = read_keyword_deck(deck_path)
    mesh = meshset(mesh_input)
    fields = Fields.allocate(mesh)
    state = RunState.from_config(cfg)

    # Initialization
    setup_initial_state(mesh, fields, cfg, state)
    apply_boundary_conditions(mesh, fields, state,
                              iwl=cfg.iwl, iwr=cfg.iwr,
                              iwt=cfg.iwt, iwb=cfg.iwb)
    free_surface_bc(mesh, fields, state, cfg)

    print(f"Input deck: {deck_path}")
    print(f"Output directory: {outdir}")
    print(f"Mesh: ibar={mesh.ibar}, jbar={mesh.jbar}, "
          f"imax={mesh.imax}, jmax={mesh.jmax}")
    print(f"Initial: Max |U|={abs(fields.U).max():.3e}, "
          f"Max |V|={abs(fields.V).max():.3e}, "
          f"F range=[{fields.F.min():.3f}, {fields.F.max():.3f}]")

    # Compute BETA relaxation coefficients (after F init, before time loop)
    compute_beta(mesh, fields, state)

    # --- Fortran pre-cycle initialization: GO TO 20 path ---
    # VFCONV (cleanup only, iCYCLE<1 skips flux), then BC, then PETACAL
    vfconv(mesh, fields, state)
    apply_boundary_conditions(mesh, fields, state,
                              iwl=cfg.iwl, iwr=cfg.iwr,
                              iwt=cfg.iwt, iwb=cfg.iwb)
    free_surface_bc(mesh, fields, state, cfg)
    petacal(mesh, fields, state, cfg)

    # --- Plot timing: mirrors Fortran pltst1/pltdt1/pltst2/pltdt2 logic ---
    pltdt = cfg.pltdt1
    twplt = cfg.pltst1
    em6 = 1.0e-6

    # Initial VTK snapshot (step=0, same as Fortran iCYCLE<=0 branch)
    write_vtk_snapshot(outdir, mesh, fields, state, step=0, prefix=cfg.prefix)
    twplt = twplt + pltdt
    if twplt + em6 > cfg.pltst2 and pltdt == cfg.pltdt1:
        pltdt = cfg.pltdt2

    # --- Time loop ---
    # Loop order matches Fortran:
    #   copy UN=U etc → zero U,V → tilde → BC → pressit → vfconv → BC → petacal
    #   → VTK write → copy UN=U etc → deltadj → advance time
    print_every = 50
    step = 0

    # Initial copy for first iteration (Fortran does this at top of loop label 10)
    fields.UN[:] = fields.U
    fields.VN[:] = fields.V
    fields.PN[:] = fields.P
    fields.FN[:] = fields.F
    fields.U.fill(0.0)
    fields.V.fill(0.0)

    while state.t < state.twfin:
        # Provisional velocity update
        tilde_step(mesh, fields, state)

        # Two-material surface tension (Fortran: between TILDE and BC)
        if cfg.nmat == 2 and cfg.isurf10:
            tms10(mesh, fields, state)

        # Boundary conditions after tilde
        apply_boundary_conditions(mesh, fields, state,
                                  iwl=cfg.iwl, iwr=cfg.iwr,
                                  iwt=cfg.iwt, iwb=cfg.iwb)
        free_surface_bc(mesh, fields, state, cfg)

        # SOLA pressure iteration (corrects P and U,V simultaneously)
        def _bc():
            apply_boundary_conditions(mesh, fields, state,
                                      iwl=cfg.iwl, iwr=cfg.iwr,
                                      iwt=cfg.iwt, iwb=cfg.iwb)
            free_surface_bc(mesh, fields, state, cfg)
        pressure_iteration(mesh, fields, state, bc_func=_bc, bc_args=())

        # VOF advection (donor-acceptor F transport)
        vfconv(mesh, fields, state)

        # Boundary conditions after VFCONV
        apply_boundary_conditions(mesh, fields, state,
                                  iwl=cfg.iwl, iwr=cfg.iwr,
                                  iwt=cfg.iwt, iwb=cfg.iwb)
        free_surface_bc(mesh, fields, state, cfg)

        # Compute NF flags for surface cells
        petacal(mesh, fields, state, cfg)

        step += 1

        if step % print_every == 0 or state.t + state.delt >= state.twfin:
            umax = abs(fields.U).max()
            vmax = abs(fields.V).max()
            print(f"t={state.t:.4f}, step={step}, delt={state.delt:.3e}, "
                  f"Max |U|={umax:.3e}, Max |V|={vmax:.3e}, "
                  f"VCHGT={state.vchgt:.3e}")

        # VTK snapshot at plot times (mirrors Fortran logic)
        if state.t + state.delt + em6 >= twplt:
            if twplt + em6 > cfg.pltst2 and pltdt == cfg.pltdt1:
                pltdt = cfg.pltdt2
            twplt = twplt + pltdt
            write_vtk_snapshot(outdir, mesh, fields, state, step=step, prefix=cfg.prefix)

        # Copy to time-n arrays (Fortran: UN=U, VN=V, PN=P, FN=F)
        # Must happen BEFORE deltadj (which reads UN/VN for velocity derivatives)
        fields.UN[:] = fields.U
        fields.VN[:] = fields.V
        fields.PN[:] = fields.P
        fields.FN[:] = fields.F

        # Adjust time step (may rollback, halve delt, recompute BETA)
        deltadj(mesh, fields, state)

        # Advance time and cycle counter
        state.t += state.delt
        state.icyle += 1

        # Zero U,V for next tilde_step
        fields.U.fill(0.0)
        fields.V.fill(0.0)

    print("Time loop finished.")


if __name__ == "__main__":
    main()

