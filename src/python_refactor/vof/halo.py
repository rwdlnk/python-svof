# halo.py — MPI halo exchange for SOLA-VOF domain decomposition
import numpy as np


def halo_exchange(decomp, *arrays):
    """
    Exchange 1-cell-wide halos with 4 neighbors (left/right/bottom/top).
    Non-blocking MPI with explicit buffer management.

    Tags are offset by array index to avoid collisions when exchanging
    multiple arrays simultaneously.
    """
    from mpi4py import MPI
    comm = decomp.comm

    if decomp.size == 1:
        return

    PROC_NULL = MPI.PROC_NULL
    requests = []
    recv_ops = []

    for arr_idx, arr in enumerate(arrays):
        ni, nj = arr.shape
        tag_base = arr_idx * 100  # unique tag space per array

        # Left: send col 1, recv into col 0
        if decomp.left != PROC_NULL:
            sbuf = np.ascontiguousarray(arr[1, :])
            rbuf = np.empty(nj, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=decomp.left, tag=tag_base + 10))
            requests.append(comm.Irecv(rbuf, source=decomp.left, tag=tag_base + 11))
            recv_ops.append((arr, 0, None, rbuf, 'col'))

        # Right: send col ni-2, recv into col ni-1
        if decomp.right != PROC_NULL:
            sbuf = np.ascontiguousarray(arr[ni - 2, :])
            rbuf = np.empty(nj, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=decomp.right, tag=tag_base + 11))
            requests.append(comm.Irecv(rbuf, source=decomp.right, tag=tag_base + 10))
            recv_ops.append((arr, ni - 1, None, rbuf, 'col'))

        # Bottom: send row 1, recv into row 0
        if decomp.bottom != PROC_NULL:
            sbuf = np.ascontiguousarray(arr[:, 1])
            rbuf = np.empty(ni, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=decomp.bottom, tag=tag_base + 20))
            requests.append(comm.Irecv(rbuf, source=decomp.bottom, tag=tag_base + 21))
            recv_ops.append((arr, None, 0, rbuf, 'row'))

        # Top: send row nj-2, recv into row nj-1
        if decomp.top != PROC_NULL:
            sbuf = np.ascontiguousarray(arr[:, nj - 2])
            rbuf = np.empty(ni, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=decomp.top, tag=tag_base + 21))
            requests.append(comm.Irecv(rbuf, source=decomp.top, tag=tag_base + 20))
            recv_ops.append((arr, None, nj - 1, rbuf, 'row'))

    MPI.Request.Waitall(requests)

    # Unpack received buffers
    for arr, ci, cj, rbuf, kind in recv_ops:
        if kind == 'col':
            arr[ci, :] = rbuf
        else:
            arr[:, cj] = rbuf


def halo_exchange_corners(decomp, *arrays):
    """
    Full 8-neighbor halo exchange including diagonal corners.
    Needed for petacal's 3x3 F stencil.

    Strategy: first do the standard 4-neighbor exchange, then exchange
    corner values via diagonal neighbors.
    """
    from mpi4py import MPI
    comm = decomp.comm

    if decomp.size == 1:
        return

    # First do the standard 4-neighbor exchange
    halo_exchange(decomp, *arrays)

    PROC_NULL = MPI.PROC_NULL
    cx, cy = decomp.coords
    px, py = decomp.dims

    def _neighbor_rank(dx, dy):
        """Get rank of neighbor at offset (dx, dy), handling periodicity."""
        nx = cx + dx
        ny = cy + dy
        if decomp.periodic_x:
            nx = nx % px
        elif nx < 0 or nx >= px:
            return PROC_NULL
        if decomp.periodic_y:
            ny = ny % py
        elif ny < 0 or ny >= py:
            return PROC_NULL
        return comm.Get_cart_rank([nx, ny])

    diag_bl = _neighbor_rank(-1, -1)
    diag_br = _neighbor_rank(+1, -1)
    diag_tl = _neighbor_rank(-1, +1)
    diag_tr = _neighbor_rank(+1, +1)

    requests = []
    recv_ops = []

    for arr_idx, arr in enumerate(arrays):
        ni, nj = arr.shape
        tag_base = arr_idx * 100 + 1000  # separate tag space from halo_exchange

        # Bottom-left corner: send arr[1,1], recv into arr[0,0]
        if diag_bl != PROC_NULL:
            sbuf = np.array([arr[1, 1]], dtype=arr.dtype)
            rbuf = np.empty(1, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=diag_bl, tag=tag_base + 30))
            requests.append(comm.Irecv(rbuf, source=diag_bl, tag=tag_base + 33))
            recv_ops.append((arr, 0, 0, rbuf))

        # Bottom-right corner: send arr[ni-2,1], recv into arr[ni-1,0]
        if diag_br != PROC_NULL:
            sbuf = np.array([arr[ni - 2, 1]], dtype=arr.dtype)
            rbuf = np.empty(1, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=diag_br, tag=tag_base + 31))
            requests.append(comm.Irecv(rbuf, source=diag_br, tag=tag_base + 32))
            recv_ops.append((arr, ni - 1, 0, rbuf))

        # Top-left corner: send arr[1,nj-2], recv into arr[0,nj-1]
        if diag_tl != PROC_NULL:
            sbuf = np.array([arr[1, nj - 2]], dtype=arr.dtype)
            rbuf = np.empty(1, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=diag_tl, tag=tag_base + 32))
            requests.append(comm.Irecv(rbuf, source=diag_tl, tag=tag_base + 31))
            recv_ops.append((arr, 0, nj - 1, rbuf))

        # Top-right corner: send arr[ni-2,nj-2], recv into arr[ni-1,nj-1]
        if diag_tr != PROC_NULL:
            sbuf = np.array([arr[ni - 2, nj - 2]], dtype=arr.dtype)
            rbuf = np.empty(1, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=diag_tr, tag=tag_base + 33))
            requests.append(comm.Irecv(rbuf, source=diag_tr, tag=tag_base + 30))
            recv_ops.append((arr, ni - 1, nj - 1, rbuf))

    MPI.Request.Waitall(requests)

    for arr, ci, cj, rbuf in recv_ops:
        arr[ci, cj] = rbuf[0]


def halo_accumulate(decomp, *arrays):
    """
    Reverse halo exchange: send halo cell values back to the owning rank
    and ADD them to the interior cells, then zero the halos.

    Needed after vfconv, which writes flux into neighbor halos
    (e.g., F[ia, j] where ia can be in the ghost region).
    """
    from mpi4py import MPI
    comm = decomp.comm

    if decomp.size == 1:
        return

    PROC_NULL = MPI.PROC_NULL
    requests = []
    recv_ops = []

    for arr_idx, arr in enumerate(arrays):
        ni, nj = arr.shape
        tag_base = arr_idx * 100 + 2000

        # Left ghost (col 0) -> send to left neighbor, they add to their col ni-2
        # Receive from left neighbor: their right ghost -> add to our col 1
        if decomp.left != PROC_NULL:
            sbuf = np.ascontiguousarray(arr[0, :])
            rbuf = np.empty(nj, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=decomp.left, tag=tag_base + 40))
            requests.append(comm.Irecv(rbuf, source=decomp.left, tag=tag_base + 41))
            recv_ops.append((arr, 1, None, rbuf, 'col_add'))

        # Right ghost (col ni-1) -> send to right neighbor, they add to their col 1
        # Receive from right neighbor: their left ghost -> add to our col ni-2
        if decomp.right != PROC_NULL:
            sbuf = np.ascontiguousarray(arr[ni - 1, :])
            rbuf = np.empty(nj, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=decomp.right, tag=tag_base + 41))
            requests.append(comm.Irecv(rbuf, source=decomp.right, tag=tag_base + 40))
            recv_ops.append((arr, ni - 2, None, rbuf, 'col_add'))

        # Bottom ghost (row 0) -> send to bottom neighbor, they add to their row nj-2
        # Receive from bottom neighbor: their top ghost -> add to our row 1
        if decomp.bottom != PROC_NULL:
            sbuf = np.ascontiguousarray(arr[:, 0])
            rbuf = np.empty(ni, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=decomp.bottom, tag=tag_base + 50))
            requests.append(comm.Irecv(rbuf, source=decomp.bottom, tag=tag_base + 51))
            recv_ops.append((arr, None, 1, rbuf, 'row_add'))

        # Top ghost (row nj-1) -> send to top neighbor, they add to their row 1
        # Receive from top neighbor: their bottom ghost -> add to our row nj-2
        if decomp.top != PROC_NULL:
            sbuf = np.ascontiguousarray(arr[:, nj - 1])
            rbuf = np.empty(ni, dtype=arr.dtype)
            requests.append(comm.Isend(sbuf, dest=decomp.top, tag=tag_base + 51))
            requests.append(comm.Irecv(rbuf, source=decomp.top, tag=tag_base + 50))
            recv_ops.append((arr, None, nj - 2, rbuf, 'row_add'))

    MPI.Request.Waitall(requests)

    # Add received values to interior boundary cells
    for arr, ci, cj, rbuf, kind in recv_ops:
        if kind == 'col_add':
            arr[ci, :] += rbuf
        else:
            arr[:, cj] += rbuf

    # Zero the ghost cells
    for arr in arrays:
        ni, nj = arr.shape
        if decomp.left != PROC_NULL:
            arr[0, :] = 0.0
        if decomp.right != PROC_NULL:
            arr[ni - 1, :] = 0.0
        if decomp.bottom != PROC_NULL:
            arr[:, 0] = 0.0
        if decomp.top != PROC_NULL:
            arr[:, nj - 1] = 0.0
