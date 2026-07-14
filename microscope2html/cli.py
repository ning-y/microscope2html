#!/usr/bin/env python3
import argparse
import glob
import os
import sys
import subprocess
import tifffile
import numpy as np
import base64
import io
from PIL import Image
import urllib.request
import zipfile
import shutil
import math
import json

# Increase PIL limit for huge images
Image.MAX_IMAGE_PIXELS = None

def get_metadata(path):
    """Extracts pixel size and stage position from OME-TIFF metadata."""
    with tifffile.TiffFile(path) as tif:
        ome_xml = tif.ome_metadata
        px_size = 1.0
        pos_x = 0.0
        pos_y = 0.0
        
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
        
        shape = tif.pages[0].shape
        return px_size, pos_x, pos_y, shape

def generate_tile_config(files, output_file="TileConfiguration.txt"):
    """Generates the TileConfiguration.txt file for ImageJ."""
    print(f"Generating {output_file} for {len(files)} tiles...")
    
    files_data = []
    for file in files:
        px_size, x_um, y_um, _ = get_metadata(file)
        if px_size == 0: px_size = 1.0
        x_px = x_um / px_size
        y_px = y_um / px_size 
        files_data.append((os.path.abspath(file), x_px, y_px))
            
    if not files_data:
        print("Error: No file data collected.")
        return False

    min_x = min(d[1] for d in files_data)
    min_y = min(d[2] for d in files_data)
    
    if not os.path.isabs(output_file):
        output_dir = os.path.dirname(files_data[0][0])
        output_path = os.path.join(output_dir, output_file)
    else:
        output_path = output_file

    print(f"Writing configuration to {output_path}...")
    with open(output_path, "w") as f:
        f.write("dim = 2\n")
        for file_path, x, y in files_data:
            filename = os.path.basename(file_path)
            f.write(f"{filename}; ; ({x - min_x:.2f}, {y - min_y:.2f})\n")
            
    return output_path

def check_white_channel(file_path):
    """Checks if the last channel is pure white."""
    try:
        with tifffile.TiffFile(file_path) as tif:
            page = tif.pages[0]
            data = page.asarray()
            if data.ndim == 3:
                if data.shape[2] < data.shape[0] and data.shape[2] < data.shape[1]:
                    channels = data.shape[2]
                    last_channel = data[:, :, channels-1]
                else:
                    channels = data.shape[0]
                    last_channel = data[channels-1, :, :]
                    
                if channels >= 4:
                    mean_val = np.mean(last_channel)
                    if mean_val > 250:
                        print(f"Detected white channel at index {channels-1} (mean {mean_val:.2f}). Will discard.")
                        return True
    except Exception as e:
        print(f"Warning: Could not check channels: {e}")
    return False

def generate_macro(output_macro, tile_config, output_image, input_dir, compute_overlap, fix_white_channel):
    """Generates the ImageJ macro."""
    overlap_param = "compute_overlap" if compute_overlap else ""
    input_dir = os.path.abspath(input_dir)
    
    macro_content = f"""
run("Grid/Collection stitching", "type=[Positions from file] order=[Defined by TileConfiguration] directory=[{input_dir}] layout_file={tile_config} fusion_method=[Linear Blending] regression_threshold=0.30 max/avg_displacement_threshold=2.50 absolute_displacement_threshold=3.50 {overlap_param} computation_parameters=[Save memory (but be slower)] image_output=[Fuse and display]");
"""
    if fix_white_channel:
        macro_content += """
run("Split Channels");
if (isOpen("C4-Fused")) { selectWindow("C4-Fused"); close(); }
run("Merge Channels...", "c1=C1-Fused c2=C2-Fused c3=C3-Fused create");
"""
    macro_content += f"""
run("RGB Color");
saveAs("Tiff", "{os.path.abspath(output_image)}");
eval("script", "System.exit(0);");
"""
    with open(output_macro, "w") as f:
        f.write(macro_content)

def download_assets():
    """Downloads OpenSeadragon assets if missing."""
    if not os.path.exists("openseadragon-bin-5.0.0"):
        print("Downloading OpenSeadragon...")
        try:
            urllib.request.urlretrieve("https://github.com/openseadragon/openseadragon/releases/download/v5.0.0/openseadragon-bin-5.0.0.zip", "osd.zip")
            with zipfile.ZipFile("osd.zip", 'r') as zip_ref:
                zip_ref.extractall(".")
            os.remove("osd.zip")
        except Exception as e:
            print(f"Error downloading OSD: {e}")
        
    if not os.path.exists("openseadragon-scalebar.js"):
        print("Downloading Scalebar plugin...")
        try:
            urllib.request.urlretrieve("https://raw.githubusercontent.com/usnistgov/OpenSeadragonScalebar/master/openseadragon-scalebar.js", "openseadragon-scalebar.js")
        except Exception as e:
            print(f"Error downloading Scalebar: {e}")

def create_tiled_html(image_path, output_html, pixel_size_um):
    """Creates a standalone HTML viewer with embedded tile pyramid."""
    print(f"Generating tile pyramid for {image_path}...")
    download_assets()
    
    # Open Image
    try:
        img = Image.open(image_path)
        if img.mode != 'RGB':
            img = img.convert('RGB')
    except Exception as e:
        print(f"Failed to open image: {e}")
        return

    width, height = img.size
    print(f"Original Size: {width}x{height}")
    
    # Calculate Max Level
    max_dim = max(width, height)
    max_level = int(math.ceil(math.log(max_dim, 2)))
    print(f"Max Level: {max_level}")
    
    tiles = {}
    tile_size = 256
    
    # Generate tiles for each level
    current_img = img
    
    # We iterate from Max Level down to 0
    for level in range(max_level, -1, -1):
        print(f"Processing Level {level} ({current_img.size[0]}x{current_img.size[1]})...")
        cols = int(math.ceil(current_img.size[0] / tile_size))
        rows = int(math.ceil(current_img.size[1] / tile_size))
        
        for col in range(cols):
            for row in range(rows):
                # Crop
                left = col * tile_size
                top = row * tile_size
                right = min(left + tile_size, current_img.size[0])
                bottom = min(top + tile_size, current_img.size[1])
                
                tile = current_img.crop((left, top, right, bottom))
                
                # Save to base64
                buffer = io.BytesIO()
                tile.save(buffer, format="JPEG", quality=75)
                b64 = base64.b64encode(buffer.getvalue()).decode("utf-8")
                
                key = f"{level}/{col}_{row}"
                tiles[key] = f"data:image/jpeg;base64,{b64}"
        
        # Prepare next level (downsample)
        if level > 0:
            new_w = max(1, current_img.size[0] // 2)
            new_h = max(1, current_img.size[1] // 2)
            current_img = current_img.resize((new_w, new_h), Image.Resampling.BILINEAR)

    print(f"Total tiles generated: {len(tiles)}")

    # Read Scripts
    try:
        with open("openseadragon-bin-5.0.0/openseadragon.min.js", "r") as f:
            osd_script = f.read()
        with open("openseadragon-scalebar.js", "r") as f:
            scalebar_script = f.read()
    except FileNotFoundError:
        print("Error: OSD scripts not found.")
        return

    # Embed Icons
    icon_prefix = "openseadragon-bin-5.0.0/images/"
    icons = {}
    if os.path.exists(icon_prefix):
        for icon_name in os.listdir(icon_prefix):
            if icon_name.endswith(".png"):
                with open(os.path.join(icon_prefix, icon_name), "rb") as f:
                    b64 = base64.b64encode(f.read()).decode("utf-8")
                    icons[icon_name] = f"data:image/png;base64,{b64}"
    
    html_content = f"""<!DOCTYPE html>
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
            zoomIn: {{ REST: icons['zoomin_rest.png'], GROUP: icons['zoomin_grouphover.png'], HOVER: icons['zoomin_hover.png'], DOWN: icons['zoomin_pressed.png'] }},
            zoomOut: {{ REST: icons['zoomout_rest.png'], GROUP: icons['zoomout_grouphover.png'], HOVER: icons['zoomout_hover.png'], DOWN: icons['zoomout_pressed.png'] }},
            home: {{ REST: icons['home_rest.png'], GROUP: icons['home_grouphover.png'], HOVER: icons['home_hover.png'], DOWN: icons['home_pressed.png'] }},
            fullpage: {{ REST: icons['fullpage_rest.png'], GROUP: icons['fullpage_grouphover.png'], HOVER: icons['fullpage_hover.png'], DOWN: icons['fullpage_pressed.png'] }},
            rotateleft: {{ REST: icons['rotateleft_rest.png'], GROUP: icons['rotateleft_grouphover.png'], HOVER: icons['rotateleft_hover.png'], DOWN: icons['rotateleft_pressed.png'] }},
            rotateright: {{ REST: icons['rotateright_rest.png'], GROUP: icons['rotateright_grouphover.png'], HOVER: icons['rotateright_hover.png'], DOWN: icons['rotateright_pressed.png'] }},
            flip: {{ REST: icons['flip_rest.png'], GROUP: icons['flip_grouphover.png'], HOVER: icons['flip_hover.png'], DOWN: icons['flip_pressed.png'] }},
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
            pixelsPerMeter: {1000000 / pixel_size_um if pixel_size_um > 0 else 1},
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
        f.write(html_content)
    print(f"Saved standalone HTML to {output_html}")

def main():
    parser = argparse.ArgumentParser(description="Stitch TIFF tiles and generate HTML viewer.")
    parser.add_argument("--input_dir", default=".", help="Directory containing TIFF tiles")
    parser.add_argument("--pattern", default="*M*d0.TIF", help="Glob pattern for tile filenames")
    parser.add_argument("--output", default="stitched_viewer.html", help="Output filename (HTML)")
    parser.add_argument("--no_overlap", action="store_true", help="Force disable overlap computation")
    parser.add_argument("--fix_white_channel", action="store_true", default=True, help="Discard 4th white channel if detected (default: True)")
    parser.add_argument("--skip_stitching", action="store_true", help="Skip stitching if output TIFF already exists")
    
    args = parser.parse_args()
    
    if not args.output.lower().endswith(".html"):
        output_tif = args.output + ".tif"
        output_html = args.output + ".html" 
    else:
        output_tif = args.output.replace(".html", ".tif")
        output_html = args.output
        
    search_path = os.path.join(args.input_dir, args.pattern)
    files = sorted(glob.glob(search_path))
    
    if not files:
        print(f"No files found matching {search_path}")
        sys.exit(1)
        
    print(f"Found {len(files)} tiles.")
    
    should_stitch = True
    if args.skip_stitching and os.path.exists(output_tif):
        print(f"Output TIFF {output_tif} exists. Skipping stitching.")
        should_stitch = False
        
    if should_stitch:
        config_filename = "TileConfiguration.txt"
        config_path = generate_tile_config(files, config_filename)
        if not config_path:
            sys.exit(1)
            
        config_file_for_macro = os.path.basename(config_path)
        
        if args.no_overlap:
            compute_overlap = False
        else:
            compute_overlap = True
            
        fix_white = args.fix_white_channel
        if fix_white:
            fix_white = check_white_channel(files[0])
            
        macro_file = "stitch_tiles.ijm"
        generate_macro(macro_file, config_file_for_macro, output_tif, args.input_dir, compute_overlap, fix_white)
        
        print("Running Fiji...")
        cmd = ["fiji", "--headless", "--console", "-macro", macro_file]
        
        try:
            subprocess.run(cmd, check=True)
            print(f"\nStitching complete. Output saved to {output_tif}")
        except subprocess.CalledProcessError as e:
            print(f"\nError running Fiji: {e}")
            sys.exit(1)
            
        for f in [macro_file, config_path]:
            if os.path.exists(f):
                os.remove(f)
        reg_file = os.path.join(os.path.dirname(config_path), "TileConfiguration.registered.txt")
        if os.path.exists(reg_file):
            os.remove(reg_file)

    if os.path.exists(output_tif):
        px_size, _, _, _ = get_metadata(files[0])
        create_tiled_html(output_tif, output_html, px_size)
    else:
        print("Error: Stitched TIFF not found.")
        sys.exit(1)

if __name__ == "__main__":
    main()
