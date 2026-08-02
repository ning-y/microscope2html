"""Overlap-based clustering of microscope tiles into ROIs."""


def cluster_tiles(tiles, overlap_margin=0.1):
    """Group tiles into ROIs by bounding-box overlap.

    Two tiles are considered adjacent (same ROI) when their
    bounding boxes overlap, with a configurable margin to handle
    small stage-coordinate imprecision.

    Args:
        tiles: List of (path, px_size, pos_x_um, pos_y_um, width_px, height_px).
        overlap_margin: Fraction of tile diagonal to expand each bounding
            box by before testing overlap (default 0.1).

    Returns:
        List of ROIs, each ROI is a list of tile tuples.
    """
    n = len(tiles)
    if n == 0:
        return []
    if n == 1:
        return [tiles]

    # Compute pixel-space bounding boxes
    px_size = tiles[0][1]
    margin_px = overlap_margin * (tiles[0][4] ** 2 + tiles[0][5] ** 2) ** 0.5

    bboxes = []
    for t in tiles:
        x0 = t[2] / px_size - margin_px
        y0 = t[3] / px_size - margin_px
        x1 = x0 + t[4] + 2 * margin_px
        y1 = y0 + t[5] + 2 * margin_px
        bboxes.append((x0, y0, x1, y1))

    # Build adjacency graph from bounding-box intersection
    adj = {i: [] for i in range(n)}
    for i in range(n):
        xi0, yi0, xi1, yi1 = bboxes[i]
        for j in range(i + 1, n):
            xj0, yj0, xj1, yj1 = bboxes[j]
            ox1 = max(xi0, xj0)
            oy1 = max(yi0, yj0)
            ox2 = min(xi1, xj1)
            oy2 = min(yi1, yj1)
            if ox1 < ox2 and oy1 < oy2:
                adj[i].append(j)
                adj[j].append(i)

    # Find connected components (ROIs) via DFS
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
            stack.extend(adj[node])
        rois.append(component)

    return rois