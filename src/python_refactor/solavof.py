# main.py — SOLA-VOF Python driver with optional MPI domain decomposition
import argparse
import csv
import time
from pathlib import Path

from vof.mesh import meshset
from vof.fields import Fields
from vof.io import read_keyword_deck, write_vtk_snapshot
from vof.state import RunState
from vof.bc import (apply_boundary_conditions, free_surface_bc,
                    apply_boundary_conditions_mpi, free_surface_bc_mpi)
from vof.initcond import setup_initial_state
from vof.solver import (tilde_step, pressure_iteration, vfconv, vfconv_phase_b,
                        petacal, compute_beta, deltadj, tms10,
                        compute_beta_mpi, deltadj_mpi)


TIMER_KEYS = ['tilde', 'halo', 'bc', 'pressure', 'vfconv',
              'petacal', 'deltadj', 'io', 'other', 'total']


def _print_timer_summary(timers, label="Totals"):
    """Print a formatted table of timer totals."""
    print(f"\n--- Timer {label} ---")
    for k in TIMER_KEYS:
        print(f"  {k:12s}: {timers[k]:10.4f} s")
    print()


def _write_timer_csv(filepath, csv_rows):
    """Write per-step timer data to CSV."""
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['step'] + TIMER_KEYS)
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"Timer CSV written: {filepath}")


def _try_mpi():
    """Try to import mpi4py and return (comm, rank, size) or None for serial."""
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        if comm.Get_size() > 1:
            return comm, comm.Get_rank(), comm.Get_size()
    except ImportError:
        pass
    return None


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
    parser.add_argument(
        "--use-cg",
        action="store_true",
        default=False,
        help="Use CG pressure solver instead of Gauss-Seidel (NMAT=2 only)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    deck_path = Path(args.input).resolve()
    outdir = Path(args.outdir).resolve()

    mpi_info = _try_mpi()
    use_mpi = mpi_info is not None

    if use_mpi:
        comm, rank, size = mpi_info
    else:
        rank, size = 0, 1

    if rank == 0:
        outdir.mkdir(parents=True, exist_ok=True)

    # All ranks read deck and build global mesh + config
    mesh_input, cfg = read_keyword_deck(deck_path)
    global_mesh = meshset(mesh_input)
    state = RunState.from_config(cfg)

    if use_mpi:
        _main_mpi(args, deck_path, outdir, global_mesh, cfg, state, comm,
                  rank, size)
    else:
        _main_serial(args, deck_path, outdir, global_mesh, cfg, state)


def _main_serial(args, deck_path, outdir, mesh, cfg, state):
    """Original serial code path (unchanged behavior)."""
    fields = Fields.allocate(mesh)

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

    compute_beta(mesh, fields, state)

    # Pre-cycle initialization
    vfconv(mesh, fields, state)
    apply_boundary_conditions(mesh, fields, state,
                              iwl=cfg.iwl, iwr=cfg.iwr,
                              iwt=cfg.iwt, iwb=cfg.iwb)
    free_surface_bc(mesh, fields, state, cfg)
    petacal(mesh, fields, state, cfg)

    pltdt = cfg.pltdt1
    twplt = cfg.pltst1
    em6 = 1.0e-6

    write_vtk_snapshot(outdir, mesh, fields, state, step=0, prefix=cfg.prefix)
    twplt = twplt + pltdt
    if twplt + em6 > cfg.pltst2 and pltdt == cfg.pltdt1:
        pltdt = cfg.pltdt2

    print_every = 50
    step = 0

    fields.UN[:] = fields.U
    fields.VN[:] = fields.V
    fields.PN[:] = fields.P
    fields.FN[:] = fields.F
    fields.U.fill(0.0)
    fields.V.fill(0.0)

    use_cg = args.use_cg and cfg.nmat == 2

    # Timer accumulators
    timers = {k: 0.0 for k in TIMER_KEYS}
    csv_rows = []

    while state.t < state.twfin:
        t_step_start = time.perf_counter()
        st = {k: 0.0 for k in TIMER_KEYS}

        # --- tilde ---
        t0 = time.perf_counter()
        tilde_step(mesh, fields, state)
        if cfg.nmat == 2 and cfg.isurf10:
            tms10(mesh, fields, state)
        st['tilde'] = time.perf_counter() - t0

        # --- bc (post-tilde) ---
        t0 = time.perf_counter()
        apply_boundary_conditions(mesh, fields, state,
                                  iwl=cfg.iwl, iwr=cfg.iwr,
                                  iwt=cfg.iwt, iwb=cfg.iwb)
        free_surface_bc(mesh, fields, state, cfg)
        st['bc'] += time.perf_counter() - t0

        # --- pressure ---
        t0 = time.perf_counter()
        if use_cg:
            from vof.cg_pressure import cg_pressure_solve
            from vof.decomp import Decomposition
            decomp_serial = Decomposition(
                comm=None, rank=0, size=1, dims=(1, 1), coords=(0, 0),
                left=-1, right=-1, bottom=-1, top=-1,
                has_left_wall=True, has_right_wall=True,
                has_bottom_wall=True, has_top_wall=True,
                gi_start=1, gi_end=mesh.ibar + 1,
                gj_start=1, gj_end=mesh.jbar + 1,
                ibar_local=mesh.ibar, jbar_local=mesh.jbar,
                periodic_x=mesh.periodic_x, periodic_y=mesh.periodic_y,
            )

            class _FakeComm:
                def Allreduce(self, sendbuf, recvbuf, op=None):
                    recvbuf[:] = sendbuf
                def allreduce(self, val, op=None):
                    return val

            decomp_serial.comm = _FakeComm()

            def _bc_cg():
                apply_boundary_conditions(mesh, fields, state,
                                          iwl=cfg.iwl, iwr=cfg.iwr,
                                          iwt=cfg.iwt, iwb=cfg.iwb)
                free_surface_bc(mesh, fields, state, cfg)

            cg_pressure_solve(decomp_serial, mesh, fields, state,
                              bc_func=_bc_cg)
        else:
            def _bc():
                apply_boundary_conditions(mesh, fields, state,
                                          iwl=cfg.iwl, iwr=cfg.iwr,
                                          iwt=cfg.iwt, iwb=cfg.iwb)
                free_surface_bc(mesh, fields, state, cfg)
            pressure_iteration(mesh, fields, state, bc_func=_bc, bc_args=())
        st['pressure'] = time.perf_counter() - t0

        # --- vfconv ---
        t0 = time.perf_counter()
        vfconv(mesh, fields, state)
        st['vfconv'] = time.perf_counter() - t0

        # --- bc (post-vfconv) ---
        t0 = time.perf_counter()
        apply_boundary_conditions(mesh, fields, state,
                                  iwl=cfg.iwl, iwr=cfg.iwr,
                                  iwt=cfg.iwt, iwb=cfg.iwb)
        free_surface_bc(mesh, fields, state, cfg)
        st['bc'] += time.perf_counter() - t0

        # --- petacal ---
        t0 = time.perf_counter()
        petacal(mesh, fields, state, cfg)
        st['petacal'] = time.perf_counter() - t0

        step += 1

        if step % print_every == 0 or state.t + state.delt >= state.twfin:
            umax = abs(fields.U).max()
            vmax = abs(fields.V).max()
            print(f"t={state.t:.4f}, step={step}, delt={state.delt:.3e}, "
                  f"Max |U|={umax:.3e}, Max |V|={vmax:.3e}, "
                  f"VCHGT={state.vchgt:.3e}")

        # --- io ---
        t0 = time.perf_counter()
        if state.t + state.delt + em6 >= twplt:
            if twplt + em6 > cfg.pltst2 and pltdt == cfg.pltdt1:
                pltdt = cfg.pltdt2
            twplt = twplt + pltdt
            write_vtk_snapshot(outdir, mesh, fields, state, step=step,
                               prefix=cfg.prefix)
        st['io'] = time.perf_counter() - t0

        # --- other (copy + deltadj + fill) ---
        t0 = time.perf_counter()
        fields.UN[:] = fields.U
        fields.VN[:] = fields.V
        fields.PN[:] = fields.P
        fields.FN[:] = fields.F
        st['other'] = time.perf_counter() - t0

        # --- deltadj ---
        t0 = time.perf_counter()
        deltadj(mesh, fields, state)
        st['deltadj'] = time.perf_counter() - t0

        # --- other (continued) ---
        t0 = time.perf_counter()
        state.t += state.delt
        state.icyle += 1
        fields.U.fill(0.0)
        fields.V.fill(0.0)
        st['other'] += time.perf_counter() - t0

        st['total'] = time.perf_counter() - t_step_start

        # Accumulate and record
        for k in TIMER_KEYS:
            timers[k] += st[k]
        csv_rows.append({'step': step, **st})

        if step % print_every == 0 or state.t >= state.twfin:
            print(f"  timers: " + ", ".join(
                f"{k}={timers[k]:.3f}" for k in TIMER_KEYS))

    _print_timer_summary(timers)
    csv_path = outdir / f"timers_{cfg.prefix}_{mesh.ibar}x{mesh.jbar}_1ranks.csv"
    _write_timer_csv(csv_path, csv_rows)
    print("Time loop finished.")


def _main_mpi(args, deck_path, outdir, global_mesh, cfg, state, comm,
              rank, size):
    """MPI parallel code path with domain decomposition."""
    from mpi4py import MPI
    from vof.decomp import (create_decomposition, build_local_mesh,
                            scatter_initial_fields, gather_global_fields)
    from vof.halo import halo_exchange, halo_exchange_corners, halo_accumulate
    from vof.io import write_vtk_snapshot_mpi
    from vof.cg_pressure import cg_pressure_solve

    # Create decomposition
    decomp = create_decomposition(comm, global_mesh)

    if rank == 0:
        print(f"MPI: {size} ranks, topology {decomp.dims[0]}x{decomp.dims[1]}")
        print(f"Input deck: {deck_path}")
        print(f"Output directory: {outdir}")
        print(f"Global mesh: ibar={global_mesh.ibar}, jbar={global_mesh.jbar}")

    # Build local mesh
    local_mesh = build_local_mesh(decomp, global_mesh)

    # Allocate local fields
    local_fields = Fields.allocate(local_mesh)

    # On rank 0, set up global fields for initialization
    if rank == 0:
        global_fields = Fields.allocate(global_mesh)
        setup_initial_state(global_mesh, global_fields, cfg, state)
        apply_boundary_conditions(global_mesh, global_fields, state,
                                  iwl=cfg.iwl, iwr=cfg.iwr,
                                  iwt=cfg.iwt, iwb=cfg.iwb)
        free_surface_bc(global_mesh, global_fields, state, cfg)
    else:
        global_fields = None

    # Scatter initial fields to all ranks
    scatter_initial_fields(decomp, global_mesh, global_fields, local_fields)

    # Halo exchange to fill ghosts from neighbors
    halo_exchange(decomp,
                  local_fields.U, local_fields.V,
                  local_fields.P, local_fields.F,
                  local_fields.BETA, local_fields.PETA,
                  local_fields.NF, local_fields.PS)

    # Apply local BCs on physical boundaries
    apply_boundary_conditions_mpi(decomp, local_mesh, local_fields, state,
                                  iwl=cfg.iwl, iwr=cfg.iwr,
                                  iwt=cfg.iwt, iwb=cfg.iwb)
    free_surface_bc_mpi(decomp, local_mesh, local_fields, state, cfg)

    if rank == 0:
        print(f"Local mesh on rank 0: ibar={local_mesh.ibar}, jbar={local_mesh.jbar}")

    # Compute BETA
    compute_beta_mpi(local_mesh, local_fields, state, decomp.comm)

    # Pre-cycle initialization (vfconv is a no-op when icycle=0,
    # so no halo_accumulate needed — just refresh ghosts)
    vfconv(local_mesh, local_fields, state)
    halo_exchange(decomp, local_fields.F)
    apply_boundary_conditions_mpi(decomp, local_mesh, local_fields, state,
                                  iwl=cfg.iwl, iwr=cfg.iwr,
                                  iwt=cfg.iwt, iwb=cfg.iwb)
    free_surface_bc_mpi(decomp, local_mesh, local_fields, state, cfg)

    halo_exchange_corners(decomp,
                          local_fields.F, local_fields.P, local_fields.BETA)
    petacal(local_mesh, local_fields, state, cfg)
    halo_exchange(decomp, local_fields.NF, local_fields.PS, local_fields.PETA)

    # Plot timing
    pltdt = cfg.pltdt1
    twplt = cfg.pltst1
    em6 = 1.0e-6

    # Allocate global_fields on rank 0 for gathering (if not already)
    if rank == 0 and global_fields is None:
        global_fields = Fields.allocate(global_mesh)

    # Initial VTK snapshot
    write_vtk_snapshot_mpi(decomp, outdir, global_mesh, local_fields,
                           global_fields, state, step=0, prefix=cfg.prefix)
    twplt = twplt + pltdt
    if twplt + em6 > cfg.pltst2 and pltdt == cfg.pltdt1:
        pltdt = cfg.pltdt2

    print_every = 50
    step = 0

    # Use CG for multi-rank (GS is sequential and won't work), or if --use-cg
    use_cg = (size > 1) or args.use_cg

    # Initial copy
    local_fields.UN[:] = local_fields.U
    local_fields.VN[:] = local_fields.V
    local_fields.PN[:] = local_fields.P
    local_fields.FN[:] = local_fields.F
    local_fields.U.fill(0.0)
    local_fields.V.fill(0.0)

    # Timer accumulators
    timers = {k: 0.0 for k in TIMER_KEYS}
    csv_rows = []

    while state.t < state.twfin:
        comm.barrier()
        t_step_start = time.perf_counter()
        st = {k: 0.0 for k in TIMER_KEYS}

        # --- tilde ---
        t0 = time.perf_counter()
        tilde_step(local_mesh, local_fields, state)
        if cfg.nmat == 2 and cfg.isurf10:
            tms10(local_mesh, local_fields, state)
        st['tilde'] = time.perf_counter() - t0

        # --- halo (post-tilde) ---
        t0 = time.perf_counter()
        halo_exchange(decomp, local_fields.U, local_fields.V)
        st['halo'] += time.perf_counter() - t0

        # --- bc (post-tilde) ---
        t0 = time.perf_counter()
        apply_boundary_conditions_mpi(decomp, local_mesh, local_fields, state,
                                      iwl=cfg.iwl, iwr=cfg.iwr,
                                      iwt=cfg.iwt, iwb=cfg.iwb)
        free_surface_bc_mpi(decomp, local_mesh, local_fields, state, cfg)
        st['bc'] += time.perf_counter() - t0

        # --- pressure ---
        t0 = time.perf_counter()
        if use_cg:
            def _bc_mpi():
                halo_exchange(decomp, local_fields.U, local_fields.V,
                              local_fields.P)
                apply_boundary_conditions_mpi(
                    decomp, local_mesh, local_fields, state,
                    iwl=cfg.iwl, iwr=cfg.iwr,
                    iwt=cfg.iwt, iwb=cfg.iwb)
                free_surface_bc_mpi(decomp, local_mesh, local_fields,
                                    state, cfg)

            cg_pressure_solve(decomp, local_mesh, local_fields, state,
                              bc_func=_bc_mpi)
        else:
            def _bc():
                apply_boundary_conditions_mpi(
                    decomp, local_mesh, local_fields, state,
                    iwl=cfg.iwl, iwr=cfg.iwr,
                    iwt=cfg.iwt, iwb=cfg.iwb)
                free_surface_bc_mpi(decomp, local_mesh, local_fields,
                                    state, cfg)
            pressure_iteration(local_mesh, local_fields, state,
                               bc_func=_bc, bc_args=())
        st['pressure'] = time.perf_counter() - t0

        # --- halo (post-pressure) ---
        t0 = time.perf_counter()
        halo_exchange(decomp, local_fields.U, local_fields.V, local_fields.P)
        st['halo'] += time.perf_counter() - t0

        # --- vfconv ---
        t0 = time.perf_counter()
        F = local_fields.F
        ni, nj = F.shape
        _saved = {}
        if decomp.left >= 0:
            _saved['L'] = F[0, :].copy()
        if decomp.right >= 0:
            _saved['R'] = F[ni - 1, :].copy()
        if decomp.bottom >= 0:
            _saved['B'] = F[:, 0].copy()
        if decomp.top >= 0:
            _saved['T'] = F[:, nj - 1].copy()

        vfconv(local_mesh, local_fields, state, skip_phase_b=True)

        if 'B' in _saved:
            F[:, 0] -= _saved['B']
        if 'T' in _saved:
            F[:, nj - 1] -= _saved['T']
        j_lo_sub = 1 if 'B' in _saved else 0
        j_hi_sub = nj - 1 if 'T' in _saved else nj
        if 'L' in _saved:
            F[0, j_lo_sub:j_hi_sub] -= _saved['L'][j_lo_sub:j_hi_sub]
        if 'R' in _saved:
            F[ni - 1, j_lo_sub:j_hi_sub] -= _saved['R'][j_lo_sub:j_hi_sub]
        st['vfconv'] += time.perf_counter() - t0

        # --- halo (vfconv accumulate + exchange) ---
        t0 = time.perf_counter()
        halo_accumulate(decomp, local_fields.F)
        halo_exchange(decomp, local_fields.F)
        st['halo'] += time.perf_counter() - t0

        # --- vfconv phase B ---
        t0 = time.perf_counter()
        vfconv_phase_b(local_mesh, local_fields, state)
        st['vfconv'] += time.perf_counter() - t0

        # --- bc (post-vfconv) ---
        t0 = time.perf_counter()
        apply_boundary_conditions_mpi(decomp, local_mesh, local_fields, state,
                                      iwl=cfg.iwl, iwr=cfg.iwr,
                                      iwt=cfg.iwt, iwb=cfg.iwb)
        free_surface_bc_mpi(decomp, local_mesh, local_fields, state, cfg)
        st['bc'] += time.perf_counter() - t0

        # --- petacal (includes corner halo) ---
        t0 = time.perf_counter()
        halo_exchange_corners(decomp,
                              local_fields.F, local_fields.P,
                              local_fields.BETA)
        st['halo'] += time.perf_counter() - t0

        t0 = time.perf_counter()
        petacal(local_mesh, local_fields, state, cfg)
        st['petacal'] = time.perf_counter() - t0

        t0 = time.perf_counter()
        halo_exchange(decomp, local_fields.NF, local_fields.PS,
                      local_fields.PETA)
        st['halo'] += time.perf_counter() - t0

        step += 1

        if step % print_every == 0 or state.t + state.delt >= state.twfin:
            local_umax = abs(local_fields.U).max()
            local_vmax = abs(local_fields.V).max()
            global_umax = decomp.comm.allreduce(local_umax, op=MPI.MAX)
            global_vmax = decomp.comm.allreduce(local_vmax, op=MPI.MAX)
            if rank == 0:
                print(f"t={state.t:.4f}, step={step}, delt={state.delt:.3e}, "
                      f"Max |U|={global_umax:.3e}, Max |V|={global_vmax:.3e}, "
                      f"VCHGT={state.vchgt:.3e}")

        # --- io ---
        t0 = time.perf_counter()
        if state.t + state.delt + em6 >= twplt:
            if twplt + em6 > cfg.pltst2 and pltdt == cfg.pltdt1:
                pltdt = cfg.pltdt2
            twplt = twplt + pltdt
            write_vtk_snapshot_mpi(decomp, outdir, global_mesh, local_fields,
                                   global_fields, state, step=step,
                                   prefix=cfg.prefix)
        st['io'] = time.perf_counter() - t0

        # --- other (copy) ---
        t0 = time.perf_counter()
        local_fields.UN[:] = local_fields.U
        local_fields.VN[:] = local_fields.V
        local_fields.PN[:] = local_fields.P
        local_fields.FN[:] = local_fields.F
        st['other'] += time.perf_counter() - t0

        # --- deltadj ---
        t0 = time.perf_counter()
        deltadj_mpi(local_mesh, local_fields, state, decomp.comm)
        st['deltadj'] = time.perf_counter() - t0

        # --- other (advance + fill) ---
        t0 = time.perf_counter()
        state.t += state.delt
        state.icyle += 1
        local_fields.U.fill(0.0)
        local_fields.V.fill(0.0)
        st['other'] += time.perf_counter() - t0

        st['total'] = time.perf_counter() - t_step_start

        # Accumulate and record
        for k in TIMER_KEYS:
            timers[k] += st[k]
        csv_rows.append({'step': step, **st})

        if step % print_every == 0 or state.t >= state.twfin:
            if rank == 0:
                print(f"  timers: " + ", ".join(
                    f"{k}={timers[k]:.3f}" for k in TIMER_KEYS))

    # Final totals: reduce with MAX across ranks to report bottleneck
    import numpy as np
    local_arr = np.array([timers[k] for k in TIMER_KEYS])
    global_arr = np.empty_like(local_arr)
    comm.Allreduce(local_arr, global_arr, op=MPI.MAX)
    max_timers = {k: global_arr[i] for i, k in enumerate(TIMER_KEYS)}

    if rank == 0:
        _print_timer_summary(max_timers, label=f"Totals (max across {size} ranks)")
        csv_path = (outdir /
                    f"timers_{cfg.prefix}_{global_mesh.ibar}x{global_mesh.jbar}"
                    f"_{size}ranks.csv")
        _write_timer_csv(csv_path, csv_rows)
        print("Time loop finished.")


if __name__ == "__main__":
    main()
