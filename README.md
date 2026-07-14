# microscope2html

A CLI tool to stitch microscope TIFF tiles (OME-TIFF) into a large image and generate a standalone, single-file HTML viewer using OpenSeadragon.

## Features

-   **Stitching**: Uses ImageJ/Fiji (via headless command line) to stitch tiles based on stage coordinates.
-   **Pyramid Generation**: Creates a Deep Zoom Image (DZI) pyramid for smooth zooming of gigapixel images.
-   **Standalone HTML**: Embeds the image tiles, OpenSeadragon library, and icons directly into a single HTML file for easy sharing.
-   **Metadata Support**: Extracts pixel size from OME-TIFF metadata to display a correct scalebar.

## Requirements

-   Python 3.12+
-   [Fiji (ImageJ)](https://fiji.sc/) installed and available on `PATH`.
-   Python dependencies: `numpy`, `tifffile`, `Pillow`

## Installation

### Via pip (any platform)

```bash
pip install .
```

Make sure Fiji is available on your `PATH`.

### Via Nix (NixOS / Nix package manager)

```bash
nix develop
python -m microscope2html.cli --help
```

This drops you into a shell with Python 3.12, all dependencies, and Fiji. See `flake.nix` for details.

## Usage

```bash
microscope2html --input_dir /path/to/tiles --pattern "*.TIF" --output result.html
```

### Arguments

-   `--input_dir`: Directory containing the TIFF tiles (default: current directory).
-   `--pattern`: Glob pattern to match tile filenames (default: `*M*d0.TIF`).
-   `--output`: Output HTML filename (default: `stitched_viewer.html`).
-   `--no_overlap`: Disable overlap computation during stitching (faster, but less accurate if coordinates are not perfect).
-   `--fix_white_channel`: Automatically discard the 4th channel if it is detected as pure white (common in some exports).
-   `--skip_stitching`: Skip the stitching step if the intermediate stitched TIFF file already exists.

## Example

```bash
microscope2html --input_dir ./data --pattern "Slide_*.tif" --output my_slide.html
```
