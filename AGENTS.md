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

### CLI expects Fiji on PATH

The CLI calls `fiji --headless --console -macro <macro>` directly. If running outside
the Nix devShell, ensure Fiji is installed and on `PATH`.

## Commit scopes

Use these scopes in Conventional Commits:

- `cli` — CLI entry point, argument parsing, output formatting
- `stitch` — Tile configuration, ImageJ macro generation, Fiji invocation
- `grid` — Grid-consistency split detection and ROI over-merge correction
- `html` — OpenSeadragon viewer, HTML template, scalebar
- `nix` — `flake.nix`, `flake.lock`, Nix-related changes
- `docs` — README, AGENTS.md, or other documentation
- `evos` — EVOS M7000 microscope data formats, metadata parsing, and reference documentation
