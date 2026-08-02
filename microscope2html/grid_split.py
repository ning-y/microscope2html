"""Grid-consistency splitting for over-merged Fiji ROIs.

When Fiji's pairwise phase correlation pulls disconnected tissue sections
together, the registered tile positions violate the expected grid spacing.
This module detects the compression and splits tiles into correct groups
using center-to-center distance filtering.
"""

import os
import re
from collections import deque


def parse_registered(reg_file):
    """Parse TileConfiguration.registered.txt into (name, x_px, y_px) tuples.

    Returns list of (filename, x_px, y_px). Tiles with (0,0) are excluded
    (these are tiles Fiji failed to register, anchored at origin).
    """
    positions = []
    with open(reg_file) as f:
        for line in f:
            m = re.search(r"\(([0-9.-]+),\s*([0-9.-]+)\)", line)
            if not m:
                continue
            x, y = float(m.group(1)), float(m.group(2))
            name = line.split(";")[0].strip()
            positions.append((name, x, y))
    return positions


def tc_expected_canvas(tiles):
    """Compute expected canvas dimensions from stage coordinates.

    Args:
        tiles: List of (filepath, px_size, x_um, y_um, width_px, height_px).

    Returns:
        (expected_width_px, expected_height_px)
    """
    if not tiles:
        return 0, 0

    _file, px_size, _x0, _y0, tile_w, tile_h = tiles[0]
    if px_size == 0:
        px_size = 1.0

    xs = [t[2] / px_size for t in tiles]
    ys = [t[3] / px_size for t in tiles]

    return max(xs) - min(xs) + tile_w, max(ys) - min(ys) + tile_h


def detect_compression(reg_positions, tiles, tile_w=2048, tile_h=1536,
                       threshold=0.85):
    """Check if Fiji-registered canvas is significantly smaller than TC-expected.

    Uses per-axis dimension ratio (min of X and Y) to catch compression
    that area alone masks (e.g., X expansion compensating for Y collapse).

    Args:
        reg_positions: List of (name, x_px, y_px) from parse_registered.
        tiles: Original tile metadata tuples.
        tile_w, tile_h: Tile dimensions in pixels.
        threshold: Minimum axis ratio below which compression is flagged.

    Returns:
        (compressed: bool, ratio: float) — ratio is the smaller axis ratio.
    """
    # Filter to active tiles only (Fiji may exclude some)
    active = [(n, x, y) for n, x, y in reg_positions if not (x == 0.0 and y == 0.0)]

    if len(active) < 2:
        return False, 1.0

    active_names = set(n for n, _, _ in active)

    # Actual canvas from registered positions
    reg_xs = [x for _, x, _ in active]
    reg_ys = [y for _, _, y in active]
    actual_w = max(reg_xs) - min(reg_xs) + tile_w
    actual_h = max(reg_ys) - min(reg_ys) + tile_h

    # Expected canvas: only for tiles that Fiji successfully registered
    active_tiles = [t for t in tiles
                    if os.path.basename(t[0]) in active_names]
    if len(active_tiles) < 2:
        return False, 1.0

    expected_w, expected_h = tc_expected_canvas(active_tiles)

    if expected_w == 0 or expected_h == 0:
        return False, 1.0

    ratio_w = actual_w / expected_w
    ratio_h = actual_h / expected_h
    ratio = min(ratio_w, ratio_h)
    return ratio < threshold, ratio


def grid_consistency_split(reg_positions, tile_w=2048, tile_h=1536,
                           x_near0_tol=200, x_near_adj_tol=100,
                           y_near0_tol=200, y_near_adj_tol=150):
    """Split registered tiles into groups using grid-consistency filtering.

    Builds the overlap graph from registered bounding boxes, then filters
    edges by whether the center-to-center distance matches expected grid
    spacing.  Pairs with distances in "dead zones" (neither near-zero
    nor near-expected-grid-spacing) are rejected — these are the
    cross-section bridges that Fiji created.

    Expected grid spacings are discovered from the bimodal distribution
    of center-to-center distances across all overlapping pairs:
      - X-axis: ~0 px (same column) or ~1900 px (adjacent column)
      - Y-axis: ~0 px (same row)    or ~1382 px (adjacent row)

    Returns:
        List of lists, each sub-list contains (name, x_px, y_px) tuples
        for one group.
    """
    active = [(n, x, y) for n, x, y in reg_positions if not (x == 0.0 and y == 0.0)]
    n = len(active)

    if n < 2:
        return [active]

    # Discover grid spacing from data
    dists = []
    for i in range(n):
        xi, yi = active[i][1], active[i][2]
        cxi, cyi = xi + tile_w / 2, yi + tile_h / 2
        for j in range(i + 1, n):
            xj, yj = active[j][1], active[j][2]
            # Check bounding-box overlap
            if (xi < xj + tile_w and xj < xi + tile_w and
                    yi < yj + tile_h and yj < yi + tile_h):
                cxj, cyj = xj + tile_w / 2, yj + tile_h / 2
                dists.append((abs(cxi - cxj), abs(cyi - cyj)))

    if not dists:
        return [active]

    # Find X spacing: should be ~0 or ~1900
    x_dists = sorted([d[0] for d in dists])
    x_adj = _find_upper_cluster(x_dists, min_val=1000)

    # Find Y spacing: should be ~0 or ~1382
    y_dists = sorted([d[1] for d in dists])
    y_adj = _find_upper_cluster(y_dists, min_val=800)

    # If we can't find clusters, use defaults
    if x_adj is None:
        x_adj = tile_w * 0.9  # assume 10% overlap
    if y_adj is None:
        y_adj = tile_h * 0.9

    # Build filtered adjacency graph
    adj = {i: set() for i in range(n)}
    for i in range(n):
        xi, yi = active[i][1], active[i][2]
        cxi, cyi = xi + tile_w / 2, yi + tile_h / 2
        for j in range(i + 1, n):
            xj, yj = active[j][1], active[j][2]
            if (xi < xj + tile_w and xj < xi + tile_w and
                    yi < yj + tile_h and yj < yi + tile_h):
                cxj, cyj = xj + tile_w / 2, yj + tile_h / 2
                dx, dy = abs(cxi - cxj), abs(cyi - cyj)

                # Valid grid spacings
                x_ok = (dx <= x_near0_tol) or (abs(dx - x_adj) <= x_near_adj_tol)
                y_ok = (dy <= y_near0_tol) or (abs(dy - y_adj) <= y_near_adj_tol)

                if x_ok and y_ok:
                    adj[i].add(j)
                    adj[j].add(i)

    # Connected components
    visited = set()
    groups = []
    for start in range(n):
        if start in visited:
            continue
        group = []
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited:
                continue
            visited.add(v)
            group.append(v)
            stack.extend(adj[v])
        groups.append(group)

    return [[active[i] for i in g] for g in groups]


def _find_upper_cluster(sorted_values, min_val=500, bin_width=100):
    """Find the dominant upper cluster in a multi-modal distribution.

    Bins values into `bin_width`-sized buckets, finds the largest cluster
    above `min_val`, and returns its mean.  This handles the case where
    cross-section bridges form an intermediate cluster (e.g., X ≈ 930)
    that must not be confused with the legitimate grid spacing (~1900).

    Returns the mean of the largest contiguous block of non-empty bins.
    """
    if not sorted_values:
        return None

    upper = [v for v in sorted_values if v >= min_val]
    if not upper:
        return None

    # Bin the values
    lo = int(min(upper) // bin_width) * bin_width
    hi = int(max(upper) // bin_width + 1) * bin_width
    bins = {}
    for v in upper:
        b = int(v // bin_width) * bin_width
        bins[b] = bins.get(b, 0) + 1

    # Find the largest contiguous block of non-empty bins
    bin_starts = sorted(bins.keys())
    best_start = None
    best_count = 0
    i = 0
    while i < len(bin_starts):
        j = i
        total = bins[bin_starts[i]]
        while j + 1 < len(bin_starts) and bin_starts[j + 1] == bin_starts[j] + bin_width:
            j += 1
            total += bins[bin_starts[j]]
        if total > best_count:
            best_count = total
            best_start = i
            best_end = j
        i = j + 1

    if best_start is None:
        return None

    # Mean of values in the winning bins
    winning_min = bin_starts[best_start]
    winning_max = bin_starts[best_end] + bin_width
    winning_vals = [v for v in upper if winning_min <= v < winning_max]
    return sum(winning_vals) / len(winning_vals) if winning_vals else None


def split_tiles_by_group(roi_tiles, groups):
    """Map grid-consistency groups back to original tile tuples.

    Args:
        roi_tiles: List of (filepath, px_size, x_um, y_um, w, h).
        groups: Output of grid_consistency_split — lists of (name, x, y).

    Returns:
        List of lists of tile tuples, one per group.
    """
    # Build name → tile lookup
    tile_by_name = {}
    for t in roi_tiles:
        basename = t[0].split("/")[-1] if "/" in t[0] else t[0]
        tile_by_name[basename] = t

    result = []
    for g in groups:
        group_tiles = []
        for name, _x, _y in g:
            if name in tile_by_name:
                group_tiles.append(tile_by_name[name])
        if group_tiles:
            result.append(group_tiles)
    return result
