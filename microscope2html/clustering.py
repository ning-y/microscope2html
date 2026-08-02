"""Distance-based clustering of microscope tiles into ROIs."""

import math


def cluster_tiles(tiles, threshold_multiplier=2.0):
    """Group tiles into ROIs by spatial proximity.

    Two tiles belong to the same ROI if their centre-to-centre distance
    is less than threshold_multiplier × tile_diagonal.

    Args:
        tiles: List of (path, px_size, pos_x_um, pos_y_um, width_px, height_px).
        threshold_multiplier: Multiplier for the tile diagonal distance.

    Returns:
        List of ROIs, each ROI is a list of tile tuples.
    """
    n = len(tiles)
    if n == 0:
        return []
    if n == 1:
        return [tiles]

    # Compute tile centres in pixel coordinates and the tile diagonal
    centres = []
    px_size = tiles[0][1]
    w0, h0 = tiles[0][4], tiles[0][5]
    tile_diagonal = math.sqrt(w0 ** 2 + h0 ** 2)
    threshold = threshold_multiplier * tile_diagonal

    for t in tiles:
        x_px = t[2] / px_size
        y_px = t[3] / px_size
        cx = x_px + t[4] / 2.0
        cy = y_px + t[5] / 2.0
        centres.append((cx, cy))

    # Build adjacency graph
    adj = {i: [] for i in range(n)}
    for i in range(n):
        xi, yi = centres[i]
        for j in range(i + 1, n):
            xj, yj = centres[j]
            dist = math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2)
            if dist < threshold:
                adj[i].append(j)
                adj[j].append(i)

    # Find connected components (ROIs)
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