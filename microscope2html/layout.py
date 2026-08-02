"""Axis-aligned band collapse for multi-ROI layout."""


def _merge_intervals(intervals):
    """Merge overlapping intervals into contiguous bands."""
    if not intervals:
        return []
    sorted_ints = sorted(intervals, key=lambda i: i[0])
    merged = [list(sorted_ints[0])]
    for start, end in sorted_ints[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(s, e) for s, e in merged]


def _find_gaps(bands):
    """Find empty gaps between merged bands."""
    gaps = []
    for i in range(len(bands) - 1):
        gap_start = bands[i][1]
        gap_end = bands[i + 1][0]
        if gap_end > gap_start:
            gaps.append((gap_start, gap_end))
    return gaps


def compress_layout(roi_bounds, margin=100):
    """Compute compressed positions for ROIs after collapsing empty bands.

    Empty horizontal and vertical bands (regions containing no ROI pixels)
    are collapsed. A configurable margin is left between ROIs that were
    previously separated by empty space.

    Args:
        roi_bounds: List of (x, y, width, height) tuples, one per ROI.
                    x,y are the top-left pixel coordinates in the original
                    (normalized) coordinate space.
        margin: Pixel gap to leave between ROIs after band collapse.

    Returns:
        List of (new_x, new_y) positions for each ROI in the compressed canvas.
    """
    n = len(roi_bounds)
    if n <= 1:
        return [(0, 0)] * n

    # Project ROIs onto X and Y axes
    x_intervals = [(x, x + w) for x, y, w, h in roi_bounds]
    y_intervals = [(y, y + h) for x, y, w, h in roi_bounds]

    x_bands = _merge_intervals(x_intervals)
    y_bands = _merge_intervals(y_intervals)

    x_gaps = _find_gaps(x_bands)
    y_gaps = _find_gaps(y_bands)

    # For each ROI, compute the total shift from collapsed gaps to its left/above
    def _shift(coord, gaps, margin):
        s = 0.0
        for gap_start, gap_end in gaps:
            if coord >= gap_end:
                s += max(0.0, (gap_end - gap_start) - margin)
        return s

    positions = []
    for x, y, w, h in roi_bounds:
        sx = _shift(x, x_gaps, margin)
        sy = _shift(y, y_gaps, margin)
        positions.append((x - sx, y - sy))

    # Ensure origin is at (0,0) — shift everything if not
    min_nx = min(p[0] for p in positions)
    min_ny = min(p[1] for p in positions)
    if min_nx > 0 or min_ny > 0:
        positions = [(nx - min_nx, ny - min_ny) for nx, ny in positions]

    return positions


def compute_canvas_bounds(roi_bounds, positions):
    """Compute the total canvas size after compression.

    Args:
        roi_bounds: List of (x, y, width, height).
        positions: List of (new_x, new_y) from compress_layout.

    Returns:
        (canvas_width, canvas_height) in pixels.
    """
    max_x = 0
    max_y = 0
    for (x, y, w, h), (nx, ny) in zip(roi_bounds, positions):
        max_x = max(max_x, nx + w)
        max_y = max(max_y, ny + h)
    return int(max_x), int(max_y)