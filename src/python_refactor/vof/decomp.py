# decomp.py — MPI domain decomposition for SOLA-VOF
from dataclasses import dataclass, field
from typing import Optional, Tuple
import numpy as np

from .mesh import Mesh


@dataclass
class Decomposition:
    """Describes how the global domain is split across MPI ranks."""
    comm: object            # MPI communicator (Cart topology)
    rank: int
    size: int
    dims: Tuple[int, int]   # (Px, Py) process grid
    coords: Tuple[int, int] # (cx, cy) this rank's position in the grid

    # Neighbor ranks (MPI_PROC_NULL = -1 for non-existent neighbors)
    left: int
    right: int
    bottom: int
    top: int

    # Whether this rank owns a physical (global) boundary
    has_left_wall: bool
    has_right_wall: bool
    has_bottom_wall: bool
    has_top_wall: bool

    # Global index ranges owned by this rank (in 0-based global arrays,
    # referring to interior cells only, i.e. Fortran 2..IM1 = Python 1..imax-2)
    gi_start: int   # first global interior i-index owned
    gi_end: int     # one past last global interior i-index owned
    gj_start: int
    gj_end: int

    # Local mesh dimensions (interior cells, no ghosts)
    ibar_local: int
    jbar_local: int

    # Periodicity
    periodic_x: bool
    periodic_y: bool


def _optimal_dims(nprocs: int, ibar: int, jbar: int) -> Tuple[int, int]:
    """Choose Px x Py that minimizes total halo perimeter."""
    best = None
    best_cost = float('inf')
    for px in range(1, nprocs + 1):
        if nprocs % px != 0:
            continue
        py = nprocs // px
        # Each subdomain is roughly (ibar/px) x (jbar/py)
        # Halo perimeter ~ 2*(jbar/py)*px + 2*(ibar/px)*py = 2*jbar + 2*ibar
        # But communication volume is proportional to interface area:
        #   x-interfaces: (py-1 internal + periodic) * (jbar/py) per rank ~ jbar*(px-1)/px...
        # Simpler: minimize total halo cells = 2*jbar*px + 2*ibar*py (total across all interfaces)
        cost = 2 * (jbar // py) * (px - 1) + 2 * (ibar // px) * (py - 1)
        if cost < best_cost:
            best_cost = cost
            best = (px, py)
    return best


def create_decomposition(comm, global_mesh: Mesh) -> 'Decomposition':
    """
    Create a Cartesian decomposition of the global domain.

    Parameters
    ----------
    comm : MPI.Comm
        MPI communicator (typically MPI.COMM_WORLD)
    global_mesh : Mesh
        The full global mesh

    Returns
    -------
    Decomposition
    """
    from mpi4py import MPI

    size = comm.Get_size()
    rank = comm.Get_rank()

    ibar = global_mesh.ibar
    jbar = global_mesh.jbar
    periodic_x = global_mesh.periodic_x
    periodic_y = global_mesh.periodic_y

    # Choose topology
    dims = _optimal_dims(size, ibar, jbar)

    # Create Cartesian communicator
    cart_comm = comm.Create_cart(
        dims=dims,
        periods=[periodic_x, periodic_y],
        reorder=True
    )
    cart_rank = cart_comm.Get_rank()
    coords = cart_comm.Get_coords(cart_rank)
    cx, cy = coords

    # Neighbor ranks (MPI_PROC_NULL for non-existent)
    left, right = cart_comm.Shift(0, 1)
    bottom, top = cart_comm.Shift(1, 1)

    # Boundary ownership
    has_left_wall = (cx == 0) and not periodic_x
    has_right_wall = (cx == dims[0] - 1) and not periodic_x
    has_bottom_wall = (cy == 0) and not periodic_y
    has_top_wall = (cy == dims[1] - 1) and not periodic_y

    # Divide interior cells among ranks
    # Interior cells in x: indices 1..ibar (0-based), total ibar cells
    gi_start, gi_end = _divide_range(ibar, dims[0], cx)
    gj_start, gj_end = _divide_range(jbar, dims[1], cy)

    # Shift to global 0-based interior indices (interior starts at index 1)
    gi_start += 1
    gi_end += 1
    gj_start += 1
    gj_end += 1

    ibar_local = gi_end - gi_start
    jbar_local = gj_end - gj_start

    return Decomposition(
        comm=cart_comm,
        rank=cart_rank,
        size=size,
        dims=dims,
        coords=(cx, cy),
        left=left, right=right,
        bottom=bottom, top=top,
        has_left_wall=has_left_wall,
        has_right_wall=has_right_wall,
        has_bottom_wall=has_bottom_wall,
        has_top_wall=has_top_wall,
        gi_start=gi_start, gi_end=gi_end,
        gj_start=gj_start, gj_end=gj_end,
        ibar_local=ibar_local, jbar_local=jbar_local,
        periodic_x=periodic_x, periodic_y=periodic_y,
    )


def _divide_range(n: int, nparts: int, part: int) -> Tuple[int, int]:
    """Divide n items into nparts, return (start, end) for partition `part`."""
    base = n // nparts
    remainder = n % nparts
    if part < remainder:
        start = part * (base + 1)
        end = start + base + 1
    else:
        start = remainder * (base + 1) + (part - remainder) * base
        end = start + base
    return start, end


def build_local_mesh(decomp: Decomposition, global_mesh: Mesh) -> Mesh:
    """
    Extract a local Mesh for this rank's subdomain from the global mesh.

    The local mesh has 1 ghost cell on each side (just like the global mesh).
    Interior cells correspond to global indices [gi_start, gi_end) x [gj_start, gj_end).

    Local layout (ibar_local interior cells + 2 ghost = imax_local total):
      index 0         = left ghost
      indices 1..ibar_local = interior cells (maps to global gi_start..gi_end-1)
      index ibar_local+1    = right ghost

    This means local im1 = ibar_local + 1 = imax_local - 1, and
    range(1, im1) covers interior cells, exactly matching the global kernel convention.
    """
    gi_s = decomp.gi_start
    gi_e = decomp.gi_end
    gj_s = decomp.gj_start
    gj_e = decomp.gj_end
    ibar_l = decomp.ibar_local
    jbar_l = decomp.jbar_local

    # Local mesh dimensions (with 1 ghost cell per side)
    imax_l = ibar_l + 2
    jmax_l = jbar_l + 2

    # Extract x-direction arrays
    # Global x array: faces at indices 0..global_mesh.imax-1
    # Interior cell i (0-based) uses face x[i] (left) and x[i+1] (right)
    # We need faces from gi_s-1 (left ghost left face) to gi_e+1 (right ghost right face)
    # That's gi_s-1 .. gi_e+1, which is ibar_l+3 faces... but we only store imax faces.
    # Actually: local cell 0 (left ghost) needs x_face[gi_s-1], x_face[gi_s]
    #           local cell 1 (first interior) needs x_face[gi_s], x_face[gi_s+1]
    #           ...
    #           local cell ibar_l (last interior) needs x_face[gi_e-1], x_face[gi_e]
    #           local cell ibar_l+1 (right ghost) needs x_face[gi_e], x_face[gi_e+1]
    # So we need faces gi_s-1 .. gi_e+1, total = ibar_l+3 = imax_l+1 faces
    # But Mesh.x has imax entries (ibar+2 for non-periodic), which are cell-face coords.
    # Actually looking at meshset: x has imax entries which are cell faces.
    # x[0]=left ghost face, x[1..ibar+1]=physical faces. For imax=ibar+2, last index=ibar+1.
    # delx[i] = x[i+1]-x[i] for i=0..imax-2, then delx[imax-1] set separately.

    # We need to build local x, delx, etc. from the global arrays.
    # The simplest approach: slice global delx for our cells, then reconstruct x.

    # Global delx indices for our cells (including ghosts):
    # local cell 0 (left ghost) = global cell gi_s-1
    # local cell k = global cell gi_s-1+k, for k=0..imax_l-1
    gm = global_mesh

    # Handle periodic wrap-around for ghost cells
    def _gidx_x(g):
        """Map global x-index with periodic wrapping."""
        if decomp.periodic_x:
            # Interior cells are 1..ibar, ghost cell at 0 and ibar+1
            if g < 1:
                return gm.ibar  # wrap to last interior
            if g > gm.ibar:
                return 1        # wrap to first interior
        return max(0, min(g, gm.imax - 1))

    def _gidx_y(g):
        """Map global y-index with periodic wrapping."""
        if decomp.periodic_y:
            if g < 1:
                return gm.jbar
            if g > gm.jbar:
                return 1
        return max(0, min(g, gm.jmax - 1))

    # Build local delx
    delx_l = np.empty(imax_l, dtype=float)
    for k in range(imax_l):
        gi = _gidx_x(gi_s - 1 + k)
        delx_l[k] = gm.delx[gi]

    dely_l = np.empty(jmax_l, dtype=float)
    for k in range(jmax_l):
        gj = _gidx_y(gj_s - 1 + k)
        dely_l[k] = gm.dely[gj]

    # Build local x (face coordinates) by cumulative sum of delx
    # We anchor x_l[1] to global x[gi_s] (the left face of the first interior cell)
    x_l = np.empty(imax_l, dtype=float)
    x_l[1] = gm.x[gi_s] if gi_s < len(gm.x) else gm.x[-1]
    x_l[0] = x_l[1] - delx_l[0]
    for k in range(2, imax_l):
        x_l[k] = x_l[k - 1] + delx_l[k - 1]

    y_l = np.empty(jmax_l, dtype=float)
    y_l[1] = gm.y[gj_s] if gj_s < len(gm.y) else gm.y[-1]
    y_l[0] = y_l[1] - dely_l[0]
    for k in range(2, jmax_l):
        y_l[k] = y_l[k - 1] + dely_l[k - 1]

    # Derived arrays
    xi_l = 0.5 * (x_l[:-1] + x_l[1:])   # cell centers, size imax_l-1
    yj_l = 0.5 * (y_l[:-1] + y_l[1:])

    rdx_l = 1.0 / delx_l
    rdy_l = 1.0 / dely_l

    rx_l = np.zeros_like(x_l)
    nz = np.abs(x_l) > 0.0
    rx_l[nz] = 1.0 / x_l[nz]

    rxi_l = np.zeros_like(xi_l)
    nz = np.abs(xi_l) > 0.0
    rxi_l[nz] = 1.0 / xi_l[nz]

    ryj_l = np.zeros_like(yj_l)
    nz = np.abs(yj_l) > 0.0
    ryj_l[nz] = 1.0 / yj_l[nz]

    im1_l = imax_l - 1
    jm1_l = jmax_l - 1
    im2_l = imax_l - 2
    jm2_l = jmax_l - 2

    return Mesh(
        x=x_l, xi=xi_l, delx=delx_l, rdx=rdx_l, rx=rx_l, rxi=rxi_l,
        y=y_l, yj=yj_l, dely=dely_l, rdy=rdy_l, ryj=ryj_l,
        ibar=ibar_l, jbar=jbar_l,
        imax=imax_l, jmax=jmax_l,
        im1=im1_l, jm1=jm1_l, im2=im2_l, jm2=jm2_l,
        nkx=gm.nkx, nky=gm.nky,
        xl=gm.xl, xc=gm.xc, dxmn=gm.dxmn,
        nxl=gm.nxl, nxr=gm.nxr,
        yl=gm.yl, yc=gm.yc, dymn=gm.dymn,
        nyl=gm.nyl, nyr=gm.nyr,
        periodic_x=decomp.periodic_x,
        periodic_y=decomp.periodic_y,
    )


def scatter_initial_fields(decomp: Decomposition, global_mesh: Mesh,
                           global_fields, local_fields) -> None:
    """
    Distribute initial field data from rank 0's global arrays to each rank's
    local arrays. Each rank gets its subdomain slice plus ghost cells filled
    from neighbors.

    Parameters
    ----------
    decomp : Decomposition
    global_mesh : Mesh (only meaningful on rank 0)
    global_fields : Fields (only meaningful on rank 0, None on others)
    local_fields : Fields (local arrays, already allocated)
    """
    from mpi4py import MPI
    comm = decomp.comm

    # Field names to scatter (2D arrays)
    names_2d = ['U', 'V', 'P', 'F', 'UN', 'VN', 'PN', 'FN',
                'BETA', 'PETA', 'DTANTH', 'PS', 'NF']

    for name in names_2d:
        if decomp.rank == 0:
            global_arr = getattr(global_fields, name)
        else:
            global_arr = None

        local_arr = getattr(local_fields, name)
        _scatter_2d(comm, decomp, global_arr, local_arr, global_mesh)

    # PR is a small 1D array, just broadcast
    if decomp.rank == 0:
        pr = global_fields.PR.copy()
    else:
        pr = np.zeros(20, dtype=float)
    comm.Bcast(pr, root=0)
    local_fields.PR[:] = pr


def _scatter_2d(comm, decomp, global_arr, local_arr, global_mesh):
    """Scatter a 2D global array to local subarrays via rank 0."""
    from mpi4py import MPI

    gi_s = decomp.gi_start
    gi_e = decomp.gi_end
    gj_s = decomp.gj_start
    gj_e = decomp.gj_end

    if decomp.rank == 0:
        # For each rank, extract the subdomain (with ghost overlap)
        for dest_rank in range(decomp.size):
            dest_coords = comm.Get_coords(dest_rank)
            dcx, dcy = dest_coords

            # Compute global index range for dest_rank
            dgi_s, dgi_e = _divide_range(global_mesh.ibar, decomp.dims[0], dcx)
            dgi_s += 1; dgi_e += 1
            dgj_s, dgj_e = _divide_range(global_mesh.jbar, decomp.dims[1], dcy)
            dgj_s += 1; dgj_e += 1

            # Extract with 1 ghost cell on each side
            i_lo = max(dgi_s - 1, 0)
            i_hi = min(dgi_e + 1, global_mesh.imax)
            j_lo = max(dgj_s - 1, 0)
            j_hi = min(dgj_e + 1, global_mesh.jmax)

            chunk = global_arr[i_lo:i_hi, j_lo:j_hi].copy()

            if dest_rank == 0:
                # Copy directly into local array
                ni = min(chunk.shape[0], local_arr.shape[0])
                nj = min(chunk.shape[1], local_arr.shape[1])
                local_arr[:ni, :nj] = chunk[:ni, :nj]
            else:
                comm.Send(chunk, dest=dest_rank, tag=100)
    else:
        # Receive our chunk
        imax_l = decomp.ibar_local + 2
        jmax_l = decomp.jbar_local + 2
        buf = np.empty((imax_l, jmax_l), dtype=local_arr.dtype)
        comm.Recv(buf, source=0, tag=100)
        local_arr[:buf.shape[0], :buf.shape[1]] = buf


def gather_global_fields(decomp: Decomposition, global_mesh: Mesh,
                         local_fields, global_fields) -> None:
    """
    Gather local field data back to rank 0's global arrays.
    Only interior cells are gathered (ghosts are redundant).

    Parameters
    ----------
    decomp : Decomposition
    global_mesh : Mesh (only meaningful on rank 0)
    local_fields : Fields (local arrays)
    global_fields : Fields (global arrays on rank 0, None on others)
    """
    from mpi4py import MPI
    comm = decomp.comm

    names_2d = ['U', 'V', 'P', 'F', 'NF', 'PS', 'BETA', 'PETA', 'DTANTH']

    for name in names_2d:
        local_arr = getattr(local_fields, name)
        # Extract interior (skip ghost cells)
        interior = local_arr[1:1 + decomp.ibar_local,
                             1:1 + decomp.jbar_local].copy()

        if decomp.rank == 0:
            global_arr = getattr(global_fields, name)
            # Place own data
            global_arr[decomp.gi_start:decomp.gi_end,
                       decomp.gj_start:decomp.gj_end] = interior

            # Receive from other ranks
            for src_rank in range(1, decomp.size):
                src_coords = comm.Get_coords(src_rank)
                scx, scy = src_coords
                sgi_s, sgi_e = _divide_range(global_mesh.ibar, decomp.dims[0], scx)
                sgi_s += 1; sgi_e += 1
                sgj_s, sgj_e = _divide_range(global_mesh.jbar, decomp.dims[1], scy)
                sgj_s += 1; sgj_e += 1

                ni = sgi_e - sgi_s
                nj = sgj_e - sgj_s
                buf = np.empty((ni, nj), dtype=global_arr.dtype)
                comm.Recv(buf, source=src_rank, tag=200)
                global_arr[sgi_s:sgi_e, sgj_s:sgj_e] = buf
        else:
            comm.Send(interior, dest=0, tag=200)
