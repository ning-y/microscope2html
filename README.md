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
- **EVOS scan discovery** — point at an EVOS M7000 scan directory and
  raw tiles are discovered, split per slide, and assembled into an
  aligned multi-channel viewer.
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

Two modes:

```bash
# Explicit files (assumes all tiles belong to one slide)
microscope2html files *.TIF [OPTIONS]

# Whole EVOS M7000 scan directory: discovers tiles, splits slides,
# stitches every captured channel, and writes one viewer per slide
microscope2html evos scan.2026-07-11-10-09-33/ [OPTIONS]
```

Every tile **must** carry OME stage coordinates (`PositionX` / `PositionY`).

### `evos` mode

Point `evos` at an EVOS M7000 output directory.  It selects the raw
acquisition tiles (never the pre-stitched tile maps unless no raw tiles
exist), groups them by slide and pass, and stitches every captured
channel using shared registration. The viewer provides individual
channel buttons and a pseudocoloured merge. Pass `--channel` to export
only one raw channel.

Output goes to `<scan-basename>_<slide>.html` in `--output-dir`
(e.g. `scan.2026-07-11-10-09-33_bottom_slide.html`), with `_p<NN>`
appended when a slide was acquired in multiple passes.  Slides that have
no raw tiles fall back to the instrument's pre-stitched tile maps with a
warning.

| Flag | Default | Description |
|---|---|---|
| `--output-dir DIR` | `.` | Directory for the output HTML files |
| `--slide NAME` | all | Only process this slide |
| `--channel CN` | all | Export only one raw channel, e.g. `4` or `d4` |
| `--dry-run` | off | List the selected tiles, sources, and outputs without stitching |
| `--margin PX` | `100` | Pixel margin between ROIs after empty‑space collapse |
| `--fix-white-channel` / `--no-fix-white-channel` | on | Discard a spurious 4th (white) channel if detected |
| `--keep-intermediates` | off | Keep the scratch directory (tile configurations, stitched ROIs) |

### `files` mode

Tiles are passed as **positional arguments** (glob patterns are expanded
by your shell).  This mode assumes a single slide and does not choose a
channel.

| Flag | Default | Description |
|---|---|---|
| `--output FILE`, `-o` | `stitched_viewer.html` | Output HTML file |
| `--margin PX` | `100` | Pixel margin between ROIs after empty‑space collapse |
| `--fix-white-channel` / `--no-fix-white-channel` | on | Discard a spurious 4th (white) channel if detected |
| `--keep-intermediates` | off | Keep the scratch directory |

### Back-end

microscope2html uses Fiji's Grid/Collection Stitching plugin.  Fiji must
be installed and on `PATH` (or `FIJI_DIR` set to the Fiji installation).
Tile-map fallback skips Fiji entirely.

### Examples

```bash
# Every slide in an EVOS scan directory, defaults
microscope2html evos scan.2026-07-11-10-09-33/

# Preview the selection without stitching
microscope2html evos --dry-run scan.2026-07-11-10-09-33/

# Only Bottom Slide, into a viewers/ directory
microscope2html evos --slide "Bottom Slide" --output-dir viewers scan.2026-07-11-10-09-33/

# Explicit files, stitched into slide.html
microscope2html files *.TIF --output slide.html
```

Note: the scratch directory (and the composite canvas, which can be
several GB) lives under the system temp directory; set `TMPDIR` if it
needs more space.

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
