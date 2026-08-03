"""Tile-to-ROI classification for microscope tiles.

Uses column-aware nearest-neighbor adjacency on stage coordinates.
No bounding-box overlap filter.  No scanprotocol dependency.
"""

import math
import os


# ═══════════════════════════════════════════════════════════════════
# Public entry point
# ═══════════════════════════════════════════════════════════════════

def cluster_tiles(tiles):
    """Group tiles into ROIs using column-aware adjacency.

    Purely spatial: uses stage coordinates and tile dimensions only.
    No scanprotocol dependency.

    Args:
        tiles: List of (path, px_size, pos_x_um, pos_y_um,
                        width_px, height_px).

    Returns:
        List of ROIs, each a list of tile tuples.
    """
    if len(tiles) <= 1:
        return [tiles] if tiles else []

    return _cluster(tiles)


# ═══════════════════════════════════════════════════════════════════
# Core algorithm
# ═══════════════════════════════════════════════════════════════════

def _cluster(tiles):
    """Column-aware spatial clustering.

    Phase 1 — Column discovery:
        Snap tile X positions to column centers.  Discover grid step
        in Y from within-column adjacent-tile distances.

    Phase 2 — Vertical adjacency (same column):
        For each column, connect tiles that are adjacent in Y within
        the discovered grid_y tolerance.

    Phase 3 — Horizontal adjacency (adjacent columns):
        For each pair of adjacent columns (nearest X neighbor),
        connect tiles at matching Y positions within a tight
        row-matching tolerance.

    Phase 4 — Connected components → ROIs.
    """
    n = len(tiles)
    tile_w, tile_h = tiles[0][4], tiles[0][5]
    px_size = tiles[0][1]

    # ── Phase 0: Build sorted position list ─────────────────────────
    indexed = list(enumerate(tiles))  # (original_idx, tile)

    # Sort by position for column grouping
    indexed.sort(key=lambda it: (it[1][2], it[1][3]))  # (X_µm, Y_µm)

    # ── Phase 1: Column discovery ──────────────────────────────────
    columns = _discover_columns(indexed, px_size, tile_w)

    # Discover grid_y from within-column Y differences
    all_dy = []
    for col_tiles in columns.values():
        for i in range(len(col_tiles) - 1):
            dy = abs(col_tiles[i][1][3] - col_tiles[i + 1][1][3])  # µm
            all_dy.append(dy)

    if len(all_dy) < 3:
        return _cluster_overlap(tiles)

    grid_y_um, zero_tol_y_um = _discover_grid_step(all_dy, tile_h * px_size)
    if grid_y_um is None:
        return _cluster_overlap(tiles)

    # ── Phase 2: Vertical adjacency within columns ───────────────────
    adj = {i: set() for i in range(n)}

    for col_tiles in columns.values():
        # col_tiles are sorted by Y (from column discovery)
        for i in range(len(col_tiles) - 1):
            idx_a, tile_a = col_tiles[i]
            idx_b, tile_b = col_tiles[i + 1]
            dy = abs(tile_a[3] - tile_b[3])
            if abs(dy - grid_y_um) <= zero_tol_y_um:
                adj[idx_a].add(idx_b)
                adj[idx_b].add(idx_a)

    # ── Phase 3: Horizontal adjacency between compatible columns ──
    # Instead of nearest-X-neighbor (which breaks with interleaved columns),
    # check ALL column pairs for Y-range overlap and connect tiles at
    # matching Y positions.

    sorted_cols = sorted(columns.keys())

    for i in range(len(sorted_cols)):
        col_a_x = sorted_cols[i]
        col_a_tiles = columns[col_a_x]
        ya_min = min(t[1][3] for t in col_a_tiles)
        ya_max = max(t[1][3] for t in col_a_tiles)

        for j in range(i + 1, len(sorted_cols)):
            col_b_x = sorted_cols[j]
            col_b_tiles = columns[col_b_x]
            yb_min = min(t[1][3] for t in col_b_tiles)
            yb_max = max(t[1][3] for t in col_b_tiles)

            # Y ranges must overlap (shared rows).
            y_overlap = min(ya_max, yb_max) - max(ya_min, yb_min)
            min_overlap = max(grid_y_um * 0.8, tile_h * px_size * 0.25)
            if y_overlap < min_overlap:
                continue

            # Connect tiles at matching Y positions
            _connect_matching_rows(
                col_a_tiles, col_b_tiles, adj,
                grid_y_um, zero_tol_y_um
            )

    # ── Phase 4: Split along large column gaps that indicate ROI
    #    boundaries.  A gap ≥ 1.3× tile width between adjacent
    #    columns means separate ROIs. ───────────────────────────────
    tile_w_um = tile_w * px_size
    gap_threshold = tile_w_um * 1.3

    sorted_cols = sorted(columns.keys())
    col_to_indices = {}
    for col_x, col_tiles in columns.items():
        col_to_indices[col_x] = {it[0] for it in col_tiles}

    for i in range(len(sorted_cols) - 1):
        gap = sorted_cols[i + 1] - sorted_cols[i]
        if gap >= gap_threshold:
            left_indices = set()
            right_indices = set()
            for j in range(i + 1):
                left_indices |= col_to_indices.get(sorted_cols[j], set())
            for j in range(i + 1, len(sorted_cols)):
                right_indices |= col_to_indices.get(sorted_cols[j], set())
            for a in left_indices:
                to_remove = [b for b in adj.get(a, set()) if b in right_indices]
                for b in to_remove:
                    adj[a].discard(b)
                    adj[b].discard(a)

    # ── Phase 5: Connected components ──────────────────────────────
    visited = set()
    rois = []
    for i in range(n):
        if i in visited:
            continue
        component = []
        stack = [i]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            component.append(tiles[node])
            stack.extend(adj[node] - visited)
        rois.append(component)

    return rois


def _discover_columns(indexed, px_size, tile_w):
    """Group tiles into columns by snapping to common X centers.

    Returns:
        dict mapping col_x_um -> list of (idx, tile) sorted by Y.
    """
    epsilon_um = px_size * 50  # 50 px tolerance for same column

    columns = {}  # col_x_um -> [(idx, tile), ...]

    # First pass: collect tiles by raw X
    for idx, tile in indexed:
        x_um = tile[2]
        y_um = tile[3]

        # Find or create column
        found = False
        for col_x in list(columns.keys()):
            if abs(x_um - col_x) <= epsilon_um:
                columns[col_x].append((idx, tile))
                found = True
                break
        if not found:
            columns[x_um] = [(idx, tile)]

    # Sort each column by Y descending (top to bottom)
    for col_x in columns:
        columns[col_x].sort(key=lambda it: -it[1][3])

    return columns


def _discover_grid_step(distances_um, tile_dim_um):
    """Discover the grid step from within-column Y differences.

    Uses histogram peak detection on the distances.

    Returns:
        (grid_step_um, tolerance_um) or (None, None) if discovery fails.
    """
    if not distances_um:
        return None, None

    # Simple mode-based approach
    # We expect distances to cluster around grid_y
    # Group into 10 µm bins
    bin_width = 10.0
    bins = {}
    for d in distances_um:
        b = int(d // bin_width) * bin_width
        bins[b] = bins.get(b, 0) + 1

    if not bins:
        return None, None

    # Find the dominant non-zero bin
    sorted_bins = sorted(bins.items(), key=lambda x: -x[1])

    # Skip zero/near-zero distances
    for bin_center, count in sorted_bins:
        if bin_center < tile_dim_um * 0.05:  # skip ~zero
            continue
        if count >= 3:
            grid_step = bin_center + bin_width / 2.0
            tolerance = bin_width * 3  # ±30 µm  (~100 px)
            return grid_step, tolerance

    return None, None


def _connect_matching_rows(col_a, col_b, adj, grid_y_um, y_tol):
    """Connect tiles in adjacent columns that share the same row.

    'Same row' means their Y positions are within y_tol of each other,
    or their Y positions are within y_tol when snapped to the grid.
    """
    # Build Y-index for column B
    ys_b = {tile[3]: (idx, tile) for idx, tile in col_b}

    for idx_a, tile_a in col_a:
        ya = tile_a[3]

        # Check exact match (within tolerance)
        for yb in ys_b:
            if abs(ya - yb) <= y_tol:
                idx_b = ys_b[yb][0]
                adj[idx_a].add(idx_b)
                adj[idx_b].add(idx_a)
                break


def _cluster_overlap(tiles):
    """Fallback: group tiles by bounding-box overlap."""
    n = len(tiles)
    if n <= 1:
        return [tiles] if n == 1 else []

    px_size = tiles[0][1]
    t_w, t_h = tiles[0][4], tiles[0][5]

    bboxes = []
    for t in tiles:
        x0 = t[2] / px_size
        y0 = t[3] / px_size
        bboxes.append((x0, y0, x0 + t_w, y0 + t_h))

    adj = {i: set() for i in range(n)}
    for i in range(n):
        xi0, yi0, xi1, yi1 = bboxes[i]
        for j in range(i + 1, n):
            xj0, yj0, xj1, yj1 = bboxes[j]
            if xi0 < xj1 and xj0 < xi1 and yi0 < yj1 and yj0 < yi1:
                adj[i].add(j)
                adj[j].add(i)

    visited = set()
    rois = []
    for i in range(n):
        if i in visited:
            continue
        component = []
        stack = [i]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            component.append(tiles[node])
            stack.extend(adj[node] - visited)
        rois.append(component)

    return rois