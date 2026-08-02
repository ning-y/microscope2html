# microscope2html

A CLI tool that stitches microscope TIFF tiles (OME-TIFF) into a single
image and generates a standalone, single-file HTML viewer using
OpenSeadragon.

Handles **discontiguous ROIs** automatically — clusters tiles by spatial
proximity, stitches each ROI independently, collapses empty space
between ROIs while preserving relative positions, and outputs everything
in one viewer.

## Features

- **Automatic ROI detection** — groups tiles into regions of interest
  based on stage-coordinate distance.
- **Overlap computation** — corrects sub-pixel stage imprecision via
  phase correlation.
- **Empty‑space collapse** — removes unused canvas between discontiguous
  ROIs while keeping their relative layout.
- **Standalone HTML** — embeds DZI tiles, OpenSeadragon, and icons in a
  single file. No server needed.
- **Scalebar** — pixel size extracted from OME-TIFF metadata.

## Requirements

- Python 3.9+
- Python packages: `numpy`, `tifffile`, `Pillow`
- Fiji (optional, for the `fiji` backend; see [Installation](#installation-fiji))

## Installation

### pip

```bash
pip install .
```

### Nix

```bash
nix develop
python -m microscope2html.cli --help
```

### Fiji {#installation-fiji}

The default overlap-computation backend (`fiji`) requires ImageJ/Fiji
installed and on `PATH`.  If Fiji is unavailable, pass `--backend
numpy-fft` to use the built-in pure‑Python backend instead.

### Container / sandbox

Set `FIJI_DIR` to point to a Fiji installation:

```bash
export FIJI_DIR=/path/to/Fiji.app   # directory containing the 'fiji' script
```

## Usage

```bash
microscope2html *.TIF [OPTIONS]
```

Tiles are passed as **positional arguments** (glob patterns are expanded
by your shell).  Every tile **must** carry OME stage coordinates
(`PositionX` / `PositionY`).

### Options

| Flag | Default | Description |
|---|---|---|
| `--output FILE` | `viewer.html` | Output HTML file |
| `--margin PX` | `100` | Pixel margin between ROIs after empty‑space collapse |
| `--cluster-threshold N` | `2.0` | Distance multiplier for ROI clustering — *N* × tile‑diagonal |
| `--backend NAME` | `fiji` | Overlap backend.  Choices: `fiji`, `numpy-fft` |
| `--fix-white-channel` | on | Discard a spurious 4th‑channel if detected (EVOS exports) |

### Back-ends

| Backend | Requires | Notes |
|---|---|---|
| `fiji` | Fiji on `PATH` | Uses Grid/Collection Stitching plugin.  Best quality; sub‑pixel alignment. |
| `numpy-fft` | nothing extra | Phase correlation via numpy FFT.  Integer‑pixel alignment; zero new dependencies. |

### Examples

```bash
# Scan all TIFFs in the current directory, stitch, output to slide.html
microscope2html *.TIF --output slide.html

# Use the numpy-FFT backend (no Fiji needed)
microscope2html *.TIF --output slide.html --backend numpy-fft

# Increase margin between ROIs and relax clustering
microscope2html *.TIF --margin 200 --cluster-threshold 3.0
```

## How it works

1. **Metadata extraction** — reads `PhysicalSizeX`, `PositionX`/`PositionY`
   from OME‑TIFF tags.
2. **ROI clustering** — connected‑components on an adjacency graph
   (tiles within `N × tile_diagonal` are neighbours).
3. **Per‑ROI stitching** — each ROI is stitched independently using the
   selected overlap backend.
4. **Layout compression** — axis‑aligned band collapse removes empty
   columns/rows between ROIs.
5. **Composite + DZI** — ROIs are pasted into a single composite image;
   a Deep Zoom Image pyramid is generated.
6. **HTML** — everything (OSD, icons, DZI tiles) is base64‑embedded into
   one standalone file.