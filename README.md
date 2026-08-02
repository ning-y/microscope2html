# microscope2html

A CLI tool that stitches microscope TIFF tiles (OME-TIFF) into a single
image and generates a standalone, single-file HTML viewer using
OpenSeadragon.

Handles **discontiguous ROIs** automatically — clusters tiles by
bounding-box overlap, stitches each ROI independently using Fiji,
collapses empty space between ROIs while preserving relative positions,
and outputs everything in one viewer.

**Detects and corrects Fiji over-merging** — when Fiji's pairwise phase
correlation incorrectly merges disconnected tissue sections, the
pipeline detects the compression via canvas-size comparison, splits the
ROI using grid-consistency filtering, and re-stitches each sub-group
correctly.

## Features

- **Automatic ROI detection** — groups tiles into regions of interest
  based on stage-coordinate bounding-box overlap.
- **Fiji stitching** — Grid/Collection Stitching with pairwise phase
  correlation for sub-pixel alignment.
- **Over-merge detection & correction** — detects when Fiji incorrectly
  merges disconnected tissue sections and splits them using
  grid-consistency filtering.  100% automatic.
- **Empty‑space collapse** — removes unused canvas between discontiguous
  ROIs while keeping their relative layout.
- **Standalone HTML** — embeds DZI tiles, OpenSeadragon, and icons in a
  single file.  No server needed.
- **Scalebar** — pixel size extracted from OME-TIFF metadata.

## Requirements

- Python 3.9+
- Python packages: `numpy`, `tifffile`, `Pillow`
- Fiji (Grid/Collection Stitching plugin; see [Installation](#installation-fiji))

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

The stitching backend (`fiji`) requires ImageJ/Fiji installed.

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
| `--overlap-margin N` | `0.1` | Fraction of tile diagonal to expand bounding boxes for ROI clustering |
| `--fix-white-channel` | on | Discard a spurious 4th‑channel if detected (EVOS exports) |

### Back-end

microscope2html uses Fiji's Grid/Collection Stitching plugin for all
stitching.  Fiji must be installed and on `PATH` (or `FIJI_DIR`
environment variable set to the directory containing the Fiji
`ImageJ-linux64` binary).

### Examples

```bash
# Scan all TIFFs in the current directory, stitch, output to slide.html
microscope2html *.TIF --output slide.html

# Increase margin between ROIs and relax clustering
microscope2html *.TIF --margin 200 --overlap-margin 0.2
```

## How it works

1. **Metadata extraction** — reads `PhysicalSizeX`, `PositionX`/`PositionY`
   from OME‑TIFF tags.
2. **ROI clustering** — connected‑components on an overlap graph
   (tiles whose stage-coordinate bounding boxes overlap are neighbours).
3. **Per‑ROI stitching** — each ROI is stitched independently using Fiji's
   Grid/Collection Stitching plugin.
4. **Over-merge detection** — after stitching, the registered canvas size
   is compared to the stage-coordinate-expected size.  If the smaller
   dimension is < 85% of expected, Fiji has incorrectly merged
   disconnected tissue sections.
5. **Grid-consistency split** — the registered overlap graph is filtered by
   center-to-center tile spacing: only edges matching the expected grid
   interval (~0 or ~1900 px in X, ~0 or ~1382 px in Y) are kept.
   Cross-section bridges (dx ≈ 930 px, dy ≈ 90 px) are rejected.
   Connected components become new sub-ROIs, each re-stitched.
6. **Layout compression** — axis‑aligned band collapse removes empty
   columns/rows between ROIs.
7. **Composite + DZI** — ROIs are pasted into a single composite image;
   a Deep Zoom Image pyramid is generated.
8. **HTML** — everything (OSD, icons, DZI tiles) is base64‑embedded into
   one standalone file.