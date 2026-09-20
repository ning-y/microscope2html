#!/usr/bin/env python3
"""Stitching, compositing, and viewer generation for microscope2html.

Both CLI modes (``files`` and ``evos``) share these stages.  All
intermediates are written to a caller-supplied scratch directory: tile
configurations, Fiji macros, per-ROI stitched images, and the composite
canvas.  The source tile directory is never written to.
"""

import base64
import io
import json
import math
import os
import shutil
import subprocess
import tempfile
import urllib.request
import zipfile
from contextlib import contextmanager

import numpy as np
import tifffile
from PIL import Image

from microscope2html.clustering import cluster_tiles
from microscope2html.layout import compress_layout, compute_canvas_bounds

Image.MAX_IMAGE_PIXELS = None

FIJI_DIR = os.environ.get("FIJI_DIR", "/tmp/fiji/Fiji.app")
FIJI_CMD = os.path.join(FIJI_DIR, "ImageJ-linux64") if os.path.isdir(FIJI_DIR) else "fiji"


@contextmanager
def scratch_dir(keep=False):
    """Create a temporary scratch directory; remove it on exit unless kept."""
    path = tempfile.mkdtemp(prefix="microscope2html-")
    try:
        yield path
    finally:
        if keep:
            print(f"Intermediates kept in {path}")
        else:
            shutil.rmtree(path, ignore_errors=True)


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


def generate_tile_config(tiles, config_path, base_dir):
    """Generate TileConfiguration.txt for a single ROI.

    Tile paths are written relative to *base_dir* (the Fiji ``directory``
    parameter).  Fiji resolves each config entry against that directory,
    so absolute paths do not work, but ``..``-relative paths do.
    """
    files_data = []
    for t in tiles:
        path, px_size, x_um, y_um, _w, _h = t
        if px_size == 0:
            px_size = 1.0
        x_px = x_um / px_size
        y_px = y_um / px_size
        rel_path = os.path.relpath(os.path.abspath(path), os.path.abspath(base_dir))
        files_data.append((rel_path, x_px, y_px))

    min_x = min(d[1] for d in files_data)
    min_y = min(d[2] for d in files_data)

    with open(config_path, "w") as f:
        f.write("dim = 2\n")
        for file_path, x, y in files_data:
            f.write(f"{file_path}; ; ({x - min_x:.2f}, {y - min_y:.2f})\n")

    return config_path


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


def _preserve_roi_configs(workdir, index):
    """Copy the generic Fiji inputs aside so every ROI's remain inspectable."""
    for name in ("TileConfiguration.txt", "TileConfiguration.registered.txt",
                 "stitch_tiles.ijm"):
        src = os.path.join(workdir, name)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(workdir, f"roi_{index}_{name}"))


def stitch_roi(roi_tiles, output_tif, workdir, fix_white=False,
               keep_intermediates=False, index=None):
    """Stitch a multi-tile ROI using Fiji."""
    if not _fiji_available():
        raise RuntimeError(
            "Fiji is not available. Install Fiji (or run inside `nix develop`) "
            "to stitch tiles."
        )
    _stitch_roi_fiji(roi_tiles, output_tif, workdir, fix_white,
                     keep_intermediates, index)


def _stitch_roi_fiji(roi_tiles, output_tif, workdir, fix_white=False,
                     keep_intermediates=False, index=None):
    """Stitch a multi-tile ROI using Fiji.  All inputs live in *workdir*."""
    input_dir = os.path.abspath(workdir)
    config_path = os.path.join(input_dir, "TileConfiguration.txt")
    generate_tile_config(roi_tiles, config_path, input_dir)

    macro_file = os.path.join(input_dir, "stitch_tiles.ijm")
    config_name = os.path.basename(config_path)
    generate_macro(macro_file, config_name, output_tif, input_dir, fix_white)

    print(f"  Running Fiji for {len(roi_tiles)} tiles...")
    cmd = [FIJI_CMD, "--headless", "--console", "-macro", macro_file]
    subprocess.run(cmd, check=True)

    if keep_intermediates and index is not None:
        _preserve_roi_configs(input_dir, index)

    # Cleanup per-ROI inputs (keep the retained copies made above)
    for name in ("stitch_tiles.ijm", "TileConfiguration.txt",
                 "TileConfiguration.registered.txt"):
        f = os.path.join(input_dir, name)
        if os.path.exists(f):
            os.remove(f)


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
            raise RuntimeError(f"Error downloading OpenSeadragon: {e}")

    if not os.path.exists("openseadragon-scalebar.js"):
        print("Downloading Scalebar plugin...")
        try:
            urllib.request.urlretrieve(
                "https://raw.githubusercontent.com/usnistgov/"
                "OpenSeadragonScalebar/master/openseadragon-scalebar.js",
                "openseadragon-scalebar.js",
            )
        except Exception as e:
            raise RuntimeError(f"Error downloading Scalebar: {e}")


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


def composite_rois(stitched_images, positions, workdir):
    """Paste all stitched ROI images into a single composite canvas.

    Returns path to the composite TIFF, or None if compositing fails.
    """
    # Determine canvas bounds
    roi_bounds = []
    for i, img_path in enumerate(stitched_images):
        with Image.open(img_path) as img:
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

    composite_path = os.path.join(workdir, "_composite.tif")
    # Use tifffile for BigTIFF support (canvas may exceed 4GB)
    composite_arr = np.array(composite)
    tifffile.imwrite(composite_path, composite_arr, bigtiff=True)
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


def build_viewer(stitched_results, tiles, output_html, margin, pixel_size_um,
                 workdir):
    """Lay out stitched ROIs, composite them, and emit the HTML viewer.

    ``stitched_results`` is a list of ``(image_path, tile_group)`` where
    each tile group's x/y coordinates define the ROI origin (raw tiles
    carry stage centers, tile maps carry top-left corners), and ``tiles``
    is the flat list of descriptors used to normalize the global origin.
    """
    global_min_x = min(t[2] / t[1] for t in tiles)
    global_min_y = min(t[3] / t[1] for t in tiles)

    roi_bounds = []
    for img_path, tile_group in stitched_results:
        with Image.open(img_path) as img:
            w, h = img.size
        roi_min_x = min(t[2] / t[1] for t in tile_group) - global_min_x
        roi_min_y = min(t[3] / t[1] for t in tile_group) - global_min_y
        roi_bounds.append((roi_min_x, roi_min_y, w, h))

    positions = compress_layout(roi_bounds, margin)

    composite_path = composite_rois([r[0] for r in stitched_results], positions,
                                    workdir)
    try:
        create_tiled_html_single(composite_path, output_html, pixel_size_um)
    finally:
        if os.path.exists(composite_path):
            os.remove(composite_path)


def process_raw_tiles(tiles, output_html, margin, fix_white_setting, workdir,
                      keep_intermediates=False):
    """Cluster raw tiles into ROIs, stitch each with Fiji, and emit HTML."""
    px_size = tiles[0][1]
    print(f"Pixel size: {px_size:.4f} µm/px")

    rois = cluster_tiles(tiles)
    print(f"Detected {len(rois)} ROI(s):")
    for i, roi in enumerate(rois):
        print(f"  ROI {i}: {len(roi)} tile(s)")

    fix_white = fix_white_setting
    if fix_white:
        fix_white = check_white_channel(tiles[0][0])

    stitched_results = []  # list of (image_path, tile_group)
    for i, roi in enumerate(rois):
        print(f"Processing ROI {i} ({len(roi)} tiles)...")
        if len(roi) == 1:
            stitched_results.append((roi[0][0], roi))
            print("  Single tile — no stitching needed.")
        else:
            output_tif = os.path.join(workdir, f"roi_{i}_stitched.tif")
            stitch_roi(roi, output_tif, workdir, fix_white,
                       keep_intermediates, index=i)
            stitched_results.append((output_tif, roi))
            print(f"  → Stitched to {output_tif}")

    build_viewer(stitched_results, tiles, output_html, margin, px_size, workdir)


def process_tm_rois(tm_paths, output_html, margin, workdir):
    """Composite pre-stitched tile maps (TM mosaics) into a viewer.

    Tile maps are already ROI-sized, so Fiji is skipped.  OME
    ``PositionX``/``PositionY`` is the *center* of each mosaic: convert it
    to a top-left origin using that mosaic's own extent, and rescale every
    mosaic to the finest pixel size among them (EVOS rescales large
    mosaics on export, so TM pixel sizes are not uniform).
    """
    metas = [get_metadata(path) for path in tm_paths]
    ref_px = min(m[0] for m in metas) or 1.0
    print(f"Using {len(tm_paths)} tile-map mosaic(s); "
          f"reference scale {ref_px:.4f} µm/px")

    stitched_results = []
    descriptors = []
    for i, (path, meta) in enumerate(zip(tm_paths, metas)):
        px, x_center, y_center, w, h = meta
        top_left_x = x_center - (w * px) / 2.0
        top_left_y = y_center - (h * px) / 2.0

        if not math.isclose(px, ref_px, rel_tol=1e-6, abs_tol=1e-9):
            with Image.open(path) as img:
                if img.mode != "RGB":
                    img = img.convert("RGB")
                new_w = max(1, int(round(w * px / ref_px)))
                new_h = max(1, int(round(h * px / ref_px)))
                img = img.resize((new_w, new_h), Image.Resampling.LANCZOS)
            scaled_path = os.path.join(workdir, f"tm_{i}_scaled.tif")
            img.save(scaled_path)
            print(f"  Rescaled {os.path.basename(path)} {w}×{h} @ {px:.4f} "
                  f"→ {new_w}×{new_h} @ {ref_px:.4f} µm/px")
            path, w, h = scaled_path, new_w, new_h

        descriptor = (path, ref_px, top_left_x, top_left_y, w, h)
        descriptors.append(descriptor)
        stitched_results.append((path, [descriptor]))

    build_viewer(stitched_results, descriptors, output_html, margin, ref_px,
                 workdir)