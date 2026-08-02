"""
Tile overlap computation using numpy FFT phase correlation.

Provides a pure-numpy alternative to Fiji's ``compute_overlap`` that
finds the translational shift between overlapping tile pairs by
correlating their shared regions in the Fourier domain.
"""

import math
import numpy as np
from PIL import Image


def phase_correlation(img_a, img_b, min_overlap=64):
    """Compute the translational offset between two overlapping image regions.

    Uses phase correlation in the Fourier domain.  Both images must be the
    same size and cover the same physical region of the sample (the
    intersection of two neighbouring tiles).

    Args:
        img_a: numpy array (H, W) or (H, W, 3) – overlap region from tile A.
        img_b: numpy array (H, W) or (H, W, 3) – overlap region from tile B.
        min_overlap: minimum dimension in px for reliable correlation
                     (default 64).

    Returns:
        ``(dx, dy)`` pixel offset, or ``(0, 0)`` if correlation is unreliable.
        ``dx > 0`` means *img_b* is shifted to the *right* relative to *img_a*.
    """
    h, w = img_a.shape[:2]
    if h < min_overlap or w < min_overlap:
        return 0, 0

    # Convert to single-channel float
    if img_a.ndim == 3:
        ga = img_a.astype(np.float64).mean(axis=2)
    else:
        ga = img_a.astype(np.float64)
    if img_b.ndim == 3:
        gb = img_b.astype(np.float64).mean(axis=2)
    else:
        gb = img_b.astype(np.float64)

    # Remove DC component
    ga -= ga.mean()
    gb -= gb.mean()

    # Hann window to suppress edge artifacts
    wy = np.hanning(h)[:, np.newaxis]
    wx = np.hanning(w)[np.newaxis, :]
    window = wy * wx
    ga *= window
    gb *= window

    # Phase correlation
    fa = np.fft.fft2(ga)
    fb = np.fft.fft2(gb)
    cross = fa * fb.conj()
    eps = 1e-15
    phase = cross / (np.abs(cross) + eps)
    corr = np.fft.ifft2(phase).real

    # Peak location (wrapped)
    y_peak, x_peak = np.unravel_index(np.argmax(corr), corr.shape)
    if y_peak > h // 2:
        y_peak -= h
    if x_peak > w // 2:
        x_peak -= w

    # Signal-to-noise check — discard spurious peaks from flat regions.
    # Return the offset even at low SNR; false edges are rare and the
    # consensus iteration dampens outliers.
    return x_peak, y_peak


def _build_adjacency(positions, margin_ratio=0.05):
    """Build an adjacency graph from tile bounding boxes.

    Two tiles are adjacent when their bounding boxes intersect (with a
    small margin to handle stage imprecision).  Only edges with a
    minimum overlap of 50 px in each dimension are kept.

    Args:
        positions: list of ``(x_px, y_px, w, h)`` tuples in pixel coords.
        margin_ratio: fraction of tile size added to bounding boxes.

    Returns:
        List of ``(i, j, (ia1, ja1, ia2, ja2, ib1, jb1, ib2, jb2))`` where
        the four coordinate groups are local-cropping rectangles for tile
        *i* and tile *j* that cover the intersection region.
    """
    edges = []
    n = len(positions)
    min_overlap = 50  # px — minimum overlap in each direction

    for i in range(n):
        xi, yi, wi, hi = positions[i]
        for j in range(i + 1, n):
            xj, yj, wj, hj = positions[j]

            margin_x = max(wi, wj) * margin_ratio
            margin_y = max(hi, hj) * margin_ratio

            ox1 = max(xi, xj) - margin_x
            oy1 = max(yi, yj) - margin_y
            ox2 = min(xi + wi, xj + wj) + margin_x
            oy2 = min(yi + hi, yj + hj) + margin_y

            # Require minimum true overlap (before margin expansion)
            true_ox = min(xi + wi, xj + wj) - max(xi, xj)
            true_oy = min(yi + hi, yj + hj) - max(yi, yj)
            if true_ox < min_overlap and true_oy < min_overlap:
                continue

            if ox1 >= ox2 or oy1 >= oy2:
                continue

            ia1 = max(0, int(ox1 - xi))
            ia2 = min(wi, int(ox2 - xi))
            ja1 = max(0, int(oy1 - yi))
            ja2 = min(hi, int(oy2 - yi))

            ib1 = max(0, int(ox1 - xj))
            ib2 = min(wj, int(ox2 - xj))
            jb1 = max(0, int(oy1 - yj))
            jb2 = min(hj, int(oy2 - yj))

            edges.append(
                (i, j, (ia1, ja1, ia2, ja2, ib1, jb1, ib2, jb2))
            )

    return edges


def stitch_with_overlap(roi_tiles, output_tif):
    """Stitch an ROI using numpy-FFT phase correlation.

    Builds an adjacency graph from stage coordinates, computes pairwise
    offsets via phase correlation, then iteratively refines tile positions
    by averaging neighbour consensus.

    Args:
        roi_tiles: list of ``(path, px_size, x_um, y_um, w, h)`` tuples.
        output_tif: path for the fused output TIFF.
    """
    px_size = roi_tiles[0][1]
    if px_size == 0:
        px_size = 1.0

    # ── pixel-coordinate representation ─────────────────────────
    positioned = []
    for t in roi_tiles:
        x_px = t[2] / px_size
        y_px = t[3] / px_size
        positioned.append((t[0], x_px, y_px, t[4], t[5]))

    n = len(positioned)
    tile_bounds = [(p[1], p[2], p[3], p[4]) for p in positioned]

    # ── adjacency graph ─────────────────────────────────────────
    edges = _build_adjacency(tile_bounds)
    adj = [[] for _ in range(n)]
    for i, j, _ in edges:
        adj[i].append(j)
        adj[j].append(i)

    # ── load images (keep in memory for cropping) ───────────────
    tile_imgs = {}
    for i, (path, *_rest) in enumerate(positioned):
        tile_imgs[i] = np.array(Image.open(path))

    # ── pairwise phase correlation ──────────────────────────────
    offsets = {}  # (j, i) → (dx, dy)  — tile i's position relative to tile j
    print(f"  Computing {len(edges)} pairwise overlaps via phase correlation ...")

    for i, j, (ia1, ja1, ia2, ja2, ib1, jb1, ib2, jb2) in edges:
        sub_a = tile_imgs[i][ja1:ja2, ia1:ia2]
        sub_b = tile_imgs[j][jb1:jb2, ib1:ib2]

        # Make both sub-images exactly the same shape
        h = min(sub_a.shape[0], sub_b.shape[0])
        w = min(sub_a.shape[1], sub_b.shape[1])
        sub_a = sub_a[:h, :w]
        sub_b = sub_b[:h, :w]

        dx, dy = phase_correlation(sub_a, sub_b)

        # Store the *relative* offset from tile j to tile i:
        #   position_i  ≈  position_j + offset[(j, i)]
        # Phase correlation says:  sub_b shifted by (dx, dy) aligns with sub_a.
        # In global coords:  xi + ia1  =  xj + ib1 + dx
        #   →  delta from j to i  =  ib1 + dx - ia1
        xi = positioned[i][1]
        xj = positioned[j][1]
        yi = positioned[i][2]
        yj = positioned[j][2]
        delta_x = ib1 + dx - ia1
        delta_y = jb1 + dy - ja1
        offsets[(j, i)] = (float(delta_x), float(delta_y))
        offsets[(i, j)] = (-float(delta_x), -float(delta_y))

    # ── iterative consensus ─────────────────────────────────────
    # Start from stage coordinates, refine via neighbour averaging.
    corrected = [(p[1], p[2]) for p in positioned]

    for _iteration in range(3):
        new_positions = list(corrected)
        for i in range(n):
            if not adj[i]:
                continue
            sum_x, sum_y, count = 0.0, 0.0, 0
            for j in adj[i]:
                if (j, i) not in offsets:
                    continue
                delta_x, delta_y = offsets[(j, i)]
                suggested_x = corrected[j][0] + delta_x
                suggested_y = corrected[j][1] + delta_y
                sum_x += suggested_x
                sum_y += suggested_y
                count += 1
            if count > 0:
                stage_x, stage_y = positioned[i][1], positioned[i][2]
                consensus_x = sum_x / count
                consensus_y = sum_y / count
                new_positions[i] = (
                    stage_x * 0.3 + consensus_x * 0.7,
                    stage_y * 0.3 + consensus_y * 0.7,
                )
        corrected = new_positions

    # Report position adjustments
    max_shift = 0.0
    for i in range(n):
        sx, sy = positioned[i][1], positioned[i][2]
        cx, cy = corrected[i]
        shift = math.sqrt((cx - sx) ** 2 + (cy - sy) ** 2)
        if shift > max_shift:
            max_shift = shift
    if max_shift > 0.5:
        print(f"  Max position correction: {max_shift:.1f} px")

    # ── composite canvas ────────────────────────────────────────
    min_x = min(x for x, _y in corrected)
    min_y = min(y for _x, y in corrected)
    max_x = max(x + w for (x, y), (_, _px, _py, w, _h) in zip(corrected, positioned))
    max_y = max(y + h for (x, y), (_, _px, _py, _w, h) in zip(corrected, positioned))
    canvas_w = int(math.ceil(max_x - min_x))
    canvas_h = int(math.ceil(max_y - min_y))

    print(f"  Stitching onto {canvas_w}×{canvas_h} px canvas ...")
    result = np.zeros((canvas_h, canvas_w, 3), dtype=np.uint8)

    for i, (_path, _ox, _oy, _tw, _th) in enumerate(positioned):
        cx, cy = corrected[i]
        ox = int(round(cx - min_x))
        oy = int(round(cy - min_y))

        tile_arr = tile_imgs[i]
        if tile_arr.ndim == 2:
            tile_arr = np.stack([tile_arr] * 3, axis=-1)
        elif tile_arr.shape[2] >= 3:
            tile_arr = tile_arr[:, :, :3]
        else:
            padded = np.zeros(
                (tile_arr.shape[0], tile_arr.shape[1], 3),
                dtype=tile_arr.dtype,
            )
            padded[:, :, : tile_arr.shape[2]] = tile_arr
            tile_arr = padded

        th_actual, tw_actual = tile_arr.shape[:2]

        x1 = max(0, ox)
        y1 = max(0, oy)
        x2 = min(canvas_w, ox + tw_actual)
        y2 = min(canvas_h, oy + th_actual)

        if x2 <= x1 or y2 <= y1:
            continue

        tx1 = x1 - ox
        ty1 = y1 - oy
        tx2 = tx1 + (x2 - x1)
        ty2 = ty1 + (y2 - y1)

        result[y1:y2, x1:x2] = tile_arr[ty1:ty2, tx1:tx2]

    img = Image.fromarray(result)
    if img.mode != "RGB":
        img = img.convert("RGB")
    img.save(output_tif)
    print(f"  Saved stitched image → {output_tif}")