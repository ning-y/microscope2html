# AGENTS.md — Project context for AI coding agents

## Project

**microscope2html** — A CLI tool that stitches microscope TIFF tiles (OME-TIFF) into a
single image and generates a standalone HTML viewer using OpenSeadragon.

## Caveats

### `pip install -e .` fails inside `nix develop`

The Nix store is read-only, so `pip install -e .` fails with `OSError: [Errno 30] Read-only file system`.
Two workarounds:

1. **Run directly** (no install needed):
   ```bash
   python -m microscope2html.cli --help
   ```
2. **Create a venv** (for editable installs):
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -e .
   ```

### CLI expects Fiji

The CLI calls Fiji for stitching. Fiji can be made available via:

1. **Nix devShell** (preferred): `nix develop` provides Fiji automatically.
2. **Manual install** (outside Nix):
   ```bash
   # Install Java
   sudo apt-get install -y default-jre
   # Download and extract Fiji
   mkdir -p /tmp/fiji && cd /tmp/fiji
   wget https://downloads.imagej.net/fiji/latest/fiji-latest-linux64-jdk.zip
   unzip fiji-latest-linux64-jdk.zip
   # The CLI expects Fiji.app/ImageJ-linux64 (newer Fiji uses fiji-linux-x64)
   ln -sf /tmp/fiji/Fiji /tmp/fiji/Fiji.app
   ln -sf fiji-linux-x64 /tmp/fiji/Fiji/ImageJ-linux64
   # Set FIJI_DIR so the CLI finds it
   export FIJI_DIR=/tmp/fiji/Fiji.app
   ```

The CLI code at `FIJI_DIR` defaults to `/tmp/fiji/Fiji.app`:
```python
FIJI_DIR = os.environ.get("FIJI_DIR", "/tmp/fiji/Fiji.app")
FIJI_CMD = os.path.join(FIJI_DIR, "ImageJ-linux64") if os.path.isdir(FIJI_DIR) else "fiji"
```

## Commit scopes

Use these scopes in Conventional Commits:

- `cli` — CLI entry point, argument parsing, output formatting
- `stitch` — Tile configuration, ImageJ macro generation, Fiji invocation
- `grid` — Grid-consistency split detection and ROI over-merge correction
- `cluster` — Tile-to-ROI classification and spatial clustering
- `html` — OpenSeadragon viewer, HTML template, scalebar
- `nix` — `flake.nix`, `flake.lock`, Nix-related changes
- `docs` — README, AGENTS.md, or other documentation
- `evos` — EVOS M7000 microscope data formats, metadata parsing, and reference documentation
