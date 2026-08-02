#!/usr/bin/env python3
"""microscope2html — Stitch microscope tiles and generate an HTML viewer.

Usage:
    microscope2html *.TIF [--output viewer.html] [--margin 100] [--overlap-margin 0.1]
"""

import argparse
import base64
import io
import json
import math
import os
import subprocess
import sys
import urllib.request
import zipfile

import numpy as np
import tifffile
from PIL import Image

from microscope2html.clustering import cluster_tiles
from microscope2html.grid_split import (
    parse_registered,
    detect_compression,
    grid_consistency_split,
    split_tiles_by_group,
)
from microscope2html.layout import compress_layout, compute_canvas_bounds

Image.MAX_IMAGE_PIXELS = None


def get_metadata(path):
    """Extract pixel size, stage position, and dimensions from OME-TIFF metadata.

    Returns (px_size, pos_x_um, pos_y_um, width_px, height_px).
    Raises ValueError if no stage PositionX/PositionY is found.
    """
    with tifffile.TiffFile(path) as tif:
        ome_xml = tif.ome_metadata
        px_size = 1.0
        pos_x = None
        pos_y = None

        if ome_xml:
            for line in ome_xml.splitlines():
                if 'PhysicalSizeX="' in line:
                    try:
                        px_size = float(line.split('PhysicalSizeX="')[1].split('"')[0])
                    except (IndexError, ValueError):
                        pass
                if 'PositionX="' in line:
                    try:
                        pos_x = float(line.split('PositionX="')[1].split('"')[0])
                    except (IndexError, ValueError):
                        pass
                if 'PositionY="' in line:
                    try:
                        pos_y = float(line.split('PositionY="')[1].split('"')[0])
                    except (IndexError, ValueError):
                        pass

        if pos_x is None or pos_y is None:
            raise ValueError(
                f"No stage coordinates (PositionX/PositionY) found in {path}"
            )

        shape = tif.pages[0].shape
        height, width = shape[0], shape[1]
        return px_size, pos_x, pos_y, width, height


def check_white_channel(file_path):
    """Check if the last channel is pure white (common in EVOS exports)."""
    try:
        with tifffile.TiffFile(file_path) as tif:
            page = tif.pages[0]
            data = page.asarray()
            if data.ndim == 3:
                if data.shape[2] < data.shape[0] and data.shape[2] < data.shape[1]:
                    channels = data.shape[2]
                    last_channel = data[:, :, channels - 1]
                else:
                    channels = data.shape[0]
                    last_channel = data[channels - 1, :, :]

                if channels >= 4:
                    mean_val = np.mean(last_channel)
                    if mean_val > 250:
                        print(
                            f"Detected white channel at index {channels - 1} "
                            f"(mean {mean_val:.2f}). Will discard."
                        )
                        return True
    except Exception as e:
        print(f"Warning: Could not check channels: {e}")
    return False


def generate_tile_config(tiles, output_file):
    """Generate TileConfiguration.txt for a single ROI."""
    files_data = []
    for t in tiles:
        path, px_size, x_um, y_um, _w, _h = t
        if px_size == 0:
            px_size = 1.0
        x_px = x_um / px_size
        y_px = y_um / px_size
        files_data.append((os.path.abspath(path), x_px, y_px))

    min_x = min(d[1] for d in files_data)
    min_y = min(d[2] for d in files_data)

    with open(output_file, "w") as f:
        f.write("dim = 2\n")
        for file_path, x, y in files_data:
            filename = os.path.basename(file_path)
            f.write(f"{filename}; ; ({x - min_x:.2f}, {y - min_y:.2f})\n")

    return output_file


def generate_macro(output_macro, tile_config, output_image, input_dir,
                   fix_white_channel):
    """Generate an ImageJ macro for stitching a single ROI."""
    input_dir = os.path.abspath(input_dir)
    output_image = os.path.abspath(output_image)

    macro = f"""
run("Grid/Collection stitching", "type=[Positions from file] \
order=[Defined by TileConfiguration] \
directory=[{input_dir}] \
layout_file={tile_config} \
fusion_method=[Linear Blending] \
regression_threshold=0.30 \
max/avg_displacement_threshold=2.50 \
absolute_displacement_threshold=3.50 \
compute_overlap \
computation_parameters=[Save memory (but be slower)] \
image_output=[Fuse and display]");
"""

    if fix_white_channel:
        macro += """
run("Split Channels");
if (isOpen("C4-Fused")) { selectWindow("C4-Fused"); close(); }
run("Merge Channels...", "c1=C1-Fused c2=C2-Fused c3=C3-Fused create");
"""

    macro += f"""
run("RGB Color");
saveAs("Tiff", "{output_image}");
eval("script", "System.exit(0);");
"""
    with open(output_macro, "w") as f:
        f.write(macro)


FIJI_DIR = os.environ.get("FIJI_DIR", "/tmp/fiji/Fiji.app")
FIJI_CMD = os.path.join(FIJI_DIR, "ImageJ-linux64") if os.path.isdir(FIJI_DIR) else "fiji"


def _fiji_available():
    """Check if Fiji is available."""
    try:
        subprocess.run(
            [FIJI_CMD, "--headless", "--version"],
            capture_output=True, timeout=10
        )
        return True
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def stitch_roi(roi_tiles, output_tif, fix_white=False):
    """Stitch a multi-tile ROI using Fiji."""
    if not _fiji_available():
        print("Error: Fiji is not available. Install Fiji to continue.",
              file=sys.stderr)
        sys.exit(1)
    _stitch_roi_fiji(roi_tiles, output_tif, fix_white)


def stitch_roi_recursive(roi_tiles, output_prefix, fix_white=False,
                         tile_w=2048, tile_h=1536):
    """Stitch an ROI, detect over-merging, and recursively split if needed.

    After Fiji stitching, compares the registered canvas size to the
    stage-coordinate-expected canvas.  If significantly smaller
    (cross-section pull), splits the tiles using grid-consistency
    filtering and re-stitches each sub-group.

    Returns:
        List of (image_path, tile_group) tuples, where tile_group is the
        list of tile metadata tuples used for this stitch.
    """
    if len(roi_tiles) <= 1:
        return [(roi_tiles[0][0], roi_tiles)]

    output_tif = f"{output_prefix}_stitched.tif"
    reg_file = _stitch_roi_fiji(roi_tiles, output_tif, fix_white,
                                keep_registered=True)

    if reg_file is None:
        print(f"  Warning: no registered positions — skipping split check")
        return [(output_tif, roi_tiles)]

    # Parse registered positions
    reg_positions = parse_registered(reg_file)
    os.remove(reg_file)

    # Detect compression
    compressed, ratio = detect_compression(reg_positions, roi_tiles,
                                           tile_w, tile_h)

    if not compressed:
        print(f"  Canvas ratio {ratio:.2f} ≥ 0.85 — no compression detected")
        return [(output_tif, roi_tiles)]

    print(f"  Canvas ratio {ratio:.2f} < 0.85 — cross-section merge detected!")
    print(f"  Running grid-consistency split...")

    # Split using grid consistency
    groups = grid_consistency_split(reg_positions, tile_w, tile_h)

    if len(groups) <= 1:
        print(f"  Grid-consistency produced only 1 group — keeping as-is")
        return [(output_tif, roi_tiles)]

    print(f"  Split into {len(groups)} groups: "
          f"{[len(g) for g in groups]}")

    # Map groups back to tile tuples
    sub_rois = split_tiles_by_group(roi_tiles, groups)

    # Remove the original stitched output (it's wrong)
    if os.path.exists(output_tif):
        os.remove(output_tif)

    # Recursively stitch each sub-group
    results = []
    for gi, sub_roi in enumerate(sub_rois):
        sub_prefix = f"{output_prefix}_g{gi}"
        results.extend(
            stitch_roi_recursive(sub_roi, sub_prefix, fix_white,
                                 tile_w, tile_h)
        )

    return results


def _stitch_roi_fiji(roi_tiles, output_tif, fix_white=False, keep_registered=False):
    """Stitch a multi-tile ROI using Fiji.

    Returns the path to TileConfiguration.registered.txt if keep_registered
    is True, otherwise None.
    """
    input_dir = os.path.dirname(os.path.abspath(roi_tiles[0][0]))
    config_path = os.path.join(input_dir, "TileConfiguration.txt")
    generate_tile_config(roi_tiles, config_path)

    macro_file = "stitch_tiles.ijm"
    config_name = os.path.basename(config_path)
    generate_macro(macro_file, config_name, output_tif, input_dir, fix_white)

    print(f"  Running Fiji for {len(roi_tiles)} tiles...")
    cmd = [FIJI_CMD, "--headless", "--console", "-macro", macro_file]
    subprocess.run(cmd, check=True)

    # Cleanup intermediate files
    for f in [macro_file, config_path]:
        if os.path.exists(f):
            os.remove(f)

    reg_file = os.path.join(input_dir, "TileConfiguration.registered.txt")
    if keep_registered and os.path.exists(reg_file):
        return reg_file
    if os.path.exists(reg_file):
        os.remove(reg_file)
    return None


def download_assets():
    """Download OpenSeadragon and scalebar plugin if missing."""
    if not os.path.exists("openseadragon-bin-5.0.0"):
        print("Downloading OpenSeadragon...")
        try:
            urllib.request.urlretrieve(
                "https://github.com/openseadragon/openseadragon/releases/"
                "download/v5.0.0/openseadragon-bin-5.0.0.zip",
                "osd.zip",
            )
            with zipfile.ZipFile("osd.zip", "r") as zip_ref:
                zip_ref.extractall(".")
            os.remove("osd.zip")
        except Exception as e:
            print(f"Error downloading OpenSeadragon: {e}")
            sys.exit(1)

    if not os.path.exists("openseadragon-scalebar.js"):
        print("Downloading Scalebar plugin...")
        try:
            urllib.request.urlretrieve(
                "https://raw.githubusercontent.com/usnistgov/"
                "OpenSeadragonScalebar/master/openseadragon-scalebar.js",
                "openseadragon-scalebar.js",
            )
        except Exception as e:
            print(f"Error downloading Scalebar: {e}")
            sys.exit(1)


def generate_dzi_tiles(image_path):
    """Generate a DZI tile pyramid from a stitched image.

    Returns (tiles_dict, width, height, max_level).
    """
    img = Image.open(image_path)
    if img.mode != "RGB":
        img = img.convert("RGB")

    width, height = img.size
    max_level = int(math.ceil(math.log(max(width, height), 2)))
    tiles = {}
    tile_size = 256
    current_img = img

    for level in range(max_level, -1, -1):
        lvl_w, lvl_h = current_img.size
        cols = int(math.ceil(lvl_w / tile_size))
        rows = int(math.ceil(lvl_h / tile_size))

        for col in range(cols):
            for row in range(rows):
                left = col * tile_size
                top = row * tile_size
                right = min(left + tile_size, lvl_w)
                bottom = min(top + tile_size, lvl_h)
                tile = current_img.crop((left, top, right, bottom))
                buf = io.BytesIO()
                tile.save(buf, format="JPEG", quality=75)
                b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
                key = f"{level}/{col}_{row}"
                tiles[key] = f"data:image/jpeg;base64,{b64}"

        if level > 0:
            new_w = max(1, lvl_w // 2)
            new_h = max(1, lvl_h // 2)
            current_img = current_img.resize((new_w, new_h), Image.Resampling.BILINEAR)

    return tiles, width, height, max_level


def composite_rois(stitched_images, positions):
    """Paste all stitched ROI images into a single composite canvas.

    Returns path to the composite TIFF, or None if compositing fails.
    """
    # Determine canvas bounds
    roi_bounds = []
    for i, img_path in enumerate(stitched_images):
        img = Image.open(img_path)
        w, h = img.size
        nx, ny = positions[i]
        roi_bounds.append((int(nx), int(ny), w, h))

    canvas_w, canvas_h = compute_canvas_bounds(
        roi_bounds, [(b[0], b[1]) for b in roi_bounds]
    )

    print(f"Compositing {len(stitched_images)} ROIs onto {canvas_w}×{canvas_h} px canvas...")

    composite = Image.new("RGB", (canvas_w, canvas_h), (0, 0, 0))

    for i, img_path in enumerate(stitched_images):
        x, y, w, h = roi_bounds[i]
        roi_img = Image.open(img_path)
        if roi_img.mode != "RGB":
            roi_img = roi_img.convert("RGB")
        composite.paste(roi_img, (x, y))

    composite_path = "_composite.tif"
    composite.save(composite_path)
    print(f"  Saved composite → {composite_path}")
    return composite_path


def create_tiled_html_single(image_path, output_html, pixel_size_um):
    """Create a standalone HTML viewer with a single embedded DZI tile pyramid."""
    print(f"Generating tile pyramid for {os.path.basename(image_path)}...")
    download_assets()

    img = Image.open(image_path)
    if img.mode != "RGB":
        img = img.convert("RGB")

    width, height = img.size
    max_level = int(math.ceil(math.log(max(width, height), 2)))
    print(f"  Size: {width}×{height}, max level: {max_level}")

    tiles = {}
    tile_size = 256
    current_img = img

    for level in range(max_level, -1, -1):
        lvl_w, lvl_h = current_img.size
        cols = int(math.ceil(lvl_w / tile_size))
        rows = int(math.ceil(lvl_h / tile_size))

        for col in range(cols):
            for row in range(rows):
                left = col * tile_size
                top = row * tile_size
                right = min(left + tile_size, lvl_w)
                bottom = min(top + tile_size, lvl_h)
                tile = current_img.crop((left, top, right, bottom))
                buf = io.BytesIO()
                tile.save(buf, format="JPEG", quality=75)
                b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
                key = f"{level}/{col}_{row}"
                tiles[key] = f"data:image/jpeg;base64,{b64}"

        if level > 0:
            new_w = max(1, lvl_w // 2)
            new_h = max(1, lvl_h // 2)
            current_img = current_img.resize((new_w, new_h), Image.Resampling.BILINEAR)

    print(f"  Total DZI tiles: {len(tiles)}")

    # Read OSD scripts
    with open("openseadragon-bin-5.0.0/openseadragon.min.js") as f:
        osd_script = f.read()
    with open("openseadragon-scalebar.js") as f:
        scalebar_script = f.read()

    # Embed icons
    icon_prefix = "openseadragon-bin-5.0.0/images/"
    icons = {}
    if os.path.exists(icon_prefix):
        for icon_name in os.listdir(icon_prefix):
            if icon_name.endswith(".png"):
                with open(os.path.join(icon_prefix, icon_name), "rb") as f:
                    b64 = base64.b64encode(f.read()).decode("utf-8")
                    icons[icon_name] = f"data:image/png;base64,{b64}"

    ppm = 1000000 / pixel_size_um if pixel_size_um > 0 else 1

    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Stitched Image Viewer</title>
    <style>
        html, body {{ width: 100%; height: 100%; margin: 0; background-color: #000; }}
        #openseadragon1 {{ width: 100%; height: 100%; }}
    </style>
</head>
<body>
    <div id="openseadragon1"></div>
    <script>
    {osd_script}
    </script>
    <script>
    {scalebar_script}
    </script>
    <script>
        var icons = {json.dumps(icons)};
        var tiles = {json.dumps(tiles)};

        var navImages = {{
            zoomIn:      {{ REST: icons['zoomin_rest.png'],      GROUP: icons['zoomin_grouphover.png'],      HOVER: icons['zoomin_hover.png'],      DOWN: icons['zoomin_pressed.png'] }},
            zoomOut:     {{ REST: icons['zoomout_rest.png'],     GROUP: icons['zoomout_grouphover.png'],     HOVER: icons['zoomout_hover.png'],     DOWN: icons['zoomout_pressed.png'] }},
            home:        {{ REST: icons['home_rest.png'],        GROUP: icons['home_grouphover.png'],        HOVER: icons['home_hover.png'],        DOWN: icons['home_pressed.png'] }},
            fullpage:    {{ REST: icons['fullpage_rest.png'],    GROUP: icons['fullpage_grouphover.png'],    HOVER: icons['fullpage_hover.png'],    DOWN: icons['fullpage_pressed.png'] }},
            rotateleft:  {{ REST: icons['rotateleft_rest.png'],  GROUP: icons['rotateleft_grouphover.png'],  HOVER: icons['rotateleft_hover.png'],  DOWN: icons['rotateleft_pressed.png'] }},
            rotateright: {{ REST: icons['rotateright_rest.png'], GROUP: icons['rotateright_grouphover.png'], HOVER: icons['rotateright_hover.png'], DOWN: icons['rotateright_pressed.png'] }},
            flip:        {{ REST: icons['flip_rest.png'],        GROUP: icons['flip_grouphover.png'],        HOVER: icons['flip_hover.png'],        DOWN: icons['flip_pressed.png'] }},
        }};

        var viewer = OpenSeadragon({{
            id: "openseadragon1",
            prefixUrl: "",
            navImages: navImages,
            tileSources: {{
                width: {width},
                height: {height},
                tileSize: {tile_size},
                tileOverlap: 0,
                minLevel: 0,
                maxLevel: {max_level},
                getTileUrl: function(level, x, y) {{
                    var key = level + "/" + x + "_" + y;
                    return tiles[key];
                }}
            }}
        }});

        viewer.scalebar({{
            pixelsPerMeter: {ppm},
            xOffset: 10,
            yOffset: 10,
            barThickness: 3,
            color: "white",
            fontColor: "white",
            backgroundColor: "rgba(0, 0, 0, 0.5)"
        }});
    </script>
</body>
</html>
"""
    with open(output_html, "w") as f:
        f.write(html)
    print(f"Saved standalone HTML to {output_html}")


def main():
    parser = argparse.ArgumentParser(
        description="Stitch microscope TIFF tiles and generate an HTML viewer."
    )
    parser.add_argument(
        "files", nargs="+",
        help="OME-TIFF tile files with stage coordinate metadata",
    )
    parser.add_argument(
        "--output", "-o", default="stitched_viewer.html",
        help="Output HTML filename (default: stitched_viewer.html)",
    )
    parser.add_argument(
        "--margin", type=int, default=100,
        help="Pixel margin between ROIs after empty band collapse (default: 100)",
    )
    parser.add_argument(
        "--overlap-margin", type=float, default=0.1,
        help="Fraction of tile diagonal to expand bounding boxes for ROI "
             "clustering (default: 0.1). Higher = more tolerant of stage "
             "coordinate imprecision.",
    )
    parser.add_argument(
        "--fix-white-channel", action="store_true", default=True,
        help="Discard 4th white channel if detected (default: True)",
    )

    args = parser.parse_args()

    # ── 1. Extract metadata from all tiles ──────────────────────────
    print(f"Reading metadata from {len(args.files)} files...")
    tiles = []
    for f in args.files:
        try:
            meta = get_metadata(f)
            tiles.append((f,) + meta)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    px_size = tiles[0][1]
    print(f"Pixel size: {px_size:.4f} µm/px")

    # ── 2. Cluster into ROIs ────────────────────────────────────────
    rois = cluster_tiles(tiles, args.overlap_margin)
    print(f"Detected {len(rois)} ROI(s):")
    for i, roi in enumerate(rois):
        print(f"  ROI {i}: {len(roi)} tile(s)")

    # ── 3. White-channel detection ──────────────────────────────────
    fix_white = args.fix_white_channel
    if fix_white:
        fix_white = check_white_channel(tiles[0][0])

    # ── 4. Stitch each ROI (with recursive split detection) ────────
    tile_w = tiles[0][4]
    tile_h = tiles[0][5]
    stitched_results = []  # list of (image_path, tile_group)
    for i, roi in enumerate(rois):
        print(f"Processing ROI {i} ({len(roi)} tiles)...")
        if len(roi) == 1:
            stitched_results.append((roi[0][0], roi))
            print(f"  Single tile — no stitching needed.")
        else:
            results = stitch_roi_recursive(
                roi, f"roi_{i}", fix_white, tile_w, tile_h
            )
            stitched_results.extend(results)
            print(f"  → {len(results)} stitched image(s)")

    stitched_images = [r[0] for r in stitched_results]

    # ── 5. Compute ROI positions in global coordinate space ─────────
    global_min_x = min(t[2] / t[1] for t in tiles)
    global_min_y = min(t[3] / t[1] for t in tiles)

    roi_bounds = []
    for img_path, tile_group in stitched_results:
        img = Image.open(img_path)
        w, h = img.size
        roi_min_x = min(t[2] / t[1] for t in tile_group) - global_min_x
        roi_min_y = min(t[3] / t[1] for t in tile_group) - global_min_y
        roi_bounds.append((roi_min_x, roi_min_y, w, h))

    # ── 6. Compress layout (collapse empty bands) ───────────────────
    positions = compress_layout(roi_bounds, args.margin)

    # ── 7. Composite all ROIs into one image, then generate HTML ────
    composite_path = composite_rois(stitched_images, positions)
    if composite_path:
        create_tiled_html_single(composite_path, args.output, px_size)
        os.remove(composite_path)

    # ── 8. Cleanup intermediate stitched TIFFs ──────────────────────
    for img_path in stitched_images:
        if img_path.startswith("roi_") and img_path.endswith(".tif"):
            if os.path.exists(img_path):
                os.remove(img_path)

    print("Done.")


if __name__ == "__main__":
    main()