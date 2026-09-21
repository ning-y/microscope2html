#!/usr/bin/env python3
"""microscope2html — Stitch microscope tiles and generate an HTML viewer.

Usage:
    microscope2html files *.TIF [OPTIONS]
    microscope2html evos SCAN_DIR/ [OPTIONS]
"""

import argparse
import math
import os
import sys
import traceback

from microscope2html import __version__
from microscope2html.evos import discover_scan, parse_channel
from microscope2html.pipeline import (
    get_metadata,
    process_raw_channels,
    process_raw_tiles,
    process_tm_rois,
    scratch_dir,
)


def _add_stitch_options(parser):
    parser.add_argument(
        "--margin", type=int, default=100,
        help="Pixel margin between ROIs after empty band collapse (default: 100)",
    )
    parser.add_argument(
        "--fix-white-channel", action=argparse.BooleanOptionalAction,
        default=True,
        help="Discard a spurious 4th (white) channel if detected "
             "(default: enabled)",
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Keep the scratch directory (tile configurations, stitched ROIs)",
    )


def build_parser():
    parser = argparse.ArgumentParser(
        prog="microscope2html",
        description="Stitch microscope TIFF tiles and generate a standalone "
                    "HTML viewer.",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}",
    )
    subparsers = parser.add_subparsers(
        dest="command", required=True, metavar="{files,evos}",
    )

    p_files = subparsers.add_parser(
        "files",
        help="stitch an explicit list of OME-TIFF tiles (assumes one slide)",
        description="Stitch an explicit list of OME-TIFF tiles.  This mode "
                    "assumes all files belong to one slide and does not "
                    "select channels; use `evos` for whole EVOS scan "
                    "directories.",
    )
    p_files.add_argument(
        "files", nargs="+",
        help="OME-TIFF tile files with stage coordinate metadata",
    )
    p_files.add_argument(
        "--output", "-o", default="stitched_viewer.html",
        help="Output HTML filename (default: stitched_viewer.html)",
    )
    _add_stitch_options(p_files)

    p_evos = subparsers.add_parser(
        "evos",
        help="discover and stitch an EVOS M7000 scan directory",
        description="Discover raw tiles in an EVOS M7000 scan directory, "
                    "split them per slide, stitch every channel, and "
                    "write one standalone HTML viewer per (slide, pass).",
    )
    p_evos.add_argument(
        "directory",
        help="EVOS scan directory, e.g. scan.2026-07-11-10-09-33/",
    )
    p_evos.add_argument(
        "--output-dir", default=".", metavar="DIR",
        help="Directory for output HTML files (default: current directory)",
    )
    p_evos.add_argument(
        "--slide", metavar="NAME",
        help="Only process this slide (exact name)",
    )
    p_evos.add_argument(
        "--channel", metavar="CN",
        help="Export only this raw channel, e.g. 4 or d4 "
             "(default: all captured channels)",
    )
    p_evos.add_argument(
        "--dry-run", action="store_true",
        help="List discovered view units and outputs without stitching",
    )
    _add_stitch_options(p_evos)

    return parser


def _validate_raw_tiles(raw_paths, label):
    """Read metadata for raw tiles and enforce one pixel grid.

    Returns ``(path, px, x, y, w, h)`` tuples ready for clustering.
    Raises ValueError naming the first offending file on mismatch.
    """
    metas = [get_metadata(path) for path in raw_paths]
    px0, _x0, _y0, w0, h0 = metas[0]
    for path, (px, _x, _y, w, h) in zip(raw_paths, metas):
        if not math.isclose(px, px0, rel_tol=1e-6, abs_tol=1e-9):
            raise ValueError(
                f"{path}: pixel size {px:g} µm/px differs from {px0:g} µm/px "
                f"in the rest of {label}"
            )
        if (w, h) != (w0, h0):
            raise ValueError(
                f"{path}: dimensions {w}×{h} differ from {w0}×{h0} "
                f"in the rest of {label}"
            )
    return [(path,) + meta for path, meta in zip(raw_paths, metas)]


def _passes_by_slide(units):
    passes = {}
    for unit in units:
        passes.setdefault(unit.slide, set()).add(unit.pas)
    return passes


def _output_name(scan_base, unit, passes_by_slide):
    suffix = ""
    if len(passes_by_slide.get(unit.slide, {unit.pas})) > 1:
        suffix = f"_p{unit.pas}"
    return f"{scan_base}_{unit.slug}{suffix}.html"


def _print_dry_run(units, directory, output_dir, problems):
    scan_base = os.path.basename(os.path.normpath(directory)) or "scan"
    passes_by_slide = _passes_by_slide(units)
    print(f"Dry run: {directory}")
    for unit in units:
        kind = ("tile maps (no stitching)" if unit.use_tm else
                "raw tiles, channels " + ", ".join(
                    f"d{channel}" for channel in unit.raw_by_channel))
        print(f"\n{unit.label}: {kind}")
        paths = unit.tm_paths if unit.use_tm else [
            path for paths in unit.raw_by_channel.values() for path in paths
        ]
        for path in paths:
            print(f"  {os.path.basename(path)}")
        if unit.raw_paths and unit.tm_paths:
            print(f"  (ignoring {len(unit.tm_paths)} tile map(s))")
        output = os.path.join(output_dir,
                              _output_name(scan_base, unit, passes_by_slide))
        print(f"  → {output}")
    for problem in problems:
        print(f"\nProblem: {problem}", file=sys.stderr)


def cmd_files(args):
    print(f"Reading metadata from {len(args.files)} files...")
    tiles = []
    for f in args.files:
        try:
            tiles.append((f,) + get_metadata(f))
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            return 1

    try:
        with scratch_dir(args.keep_intermediates) as workdir:
            process_raw_tiles(tiles, args.output, args.margin,
                              args.fix_white_channel, workdir,
                              args.keep_intermediates)
    except RuntimeError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    print("Done.")
    return 0


def cmd_evos(args):
    directory = os.path.abspath(args.directory)
    if not os.path.isdir(directory):
        print(f"Error: not a directory: {args.directory}", file=sys.stderr)
        return 1

    try:
        channel = parse_channel(args.channel) if args.channel is not None else None
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    units, problems = discover_scan(directory, slide=args.slide,
                                    channel=channel)

    valid_units = []
    for unit in units:
        try:
            if unit.raw_paths:
                unit.tiles = {
                    channel: _validate_raw_tiles(
                        paths, f"{unit.label}, channel d{channel}")
                    for channel, paths in unit.raw_by_channel.items()
                }
                reference = next(iter(unit.tiles.values()))
                reference_keys = {os.path.basename(tile[0]).rsplit("d", 1)[0]
                                  for tile in reference}
                for channel, channel_tiles in unit.tiles.items():
                    keys = {os.path.basename(tile[0]).rsplit("d", 1)[0]
                            for tile in channel_tiles}
                    if keys != reference_keys:
                        raise ValueError(
                            f"{unit.label}: channel d{channel} does not have "
                            "the same fields as the reference channel"
                        )
            else:
                for path in unit.tm_paths:
                    get_metadata(path)
        except ValueError as e:
            problems.append(f"{unit.label}: {e}")
            continue
        valid_units.append(unit)

    if args.dry_run:
        _print_dry_run(valid_units, directory, args.output_dir, problems)
        return 1 if problems else 0

    if not valid_units:
        for problem in problems:
            print(f"Error: {problem}", file=sys.stderr)
        return 1

    for problem in problems:
        print(f"Warning: {problem}", file=sys.stderr)

    try:
        os.makedirs(args.output_dir, exist_ok=True)
    except OSError as e:
        print(f"Error: cannot create output directory {args.output_dir}: {e}",
              file=sys.stderr)
        return 1

    scan_base = os.path.basename(os.path.normpath(directory)) or "scan"
    passes_by_slide = _passes_by_slide(valid_units)
    failures = 0
    for unit in valid_units:
        output_html = os.path.join(
            args.output_dir, _output_name(scan_base, unit, passes_by_slide)
        )
        print(f"\n=== {unit.label} → {output_html} ===")
        if unit.use_tm:
            print("Warning: no raw tiles for this slide; falling back to "
                  "pre-stitched tile maps from the microscope.")
        elif unit.tm_paths:
            print(f"Note: using raw tiles; ignoring "
                  f"{len(unit.tm_paths)} tile map(s).")
        try:
            with scratch_dir(args.keep_intermediates) as workdir:
                if unit.use_tm:
                    process_tm_rois(unit.tm_paths, output_html, args.margin,
                                    workdir)
                else:
                    process_raw_channels(unit.tiles, output_html, args.margin,
                                         args.fix_white_channel, workdir,
                                         args.keep_intermediates)
        except Exception as e:  # keep other units going
            failures += 1
            print(f"Error: {type(e).__name__}: {e}", file=sys.stderr)
            if not isinstance(e, RuntimeError):
                traceback.print_exc()

    if failures or problems:
        return 1
    print("\nDone.")
    return 0


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "files":
        return cmd_files(args)
    if args.command == "evos":
        return cmd_evos(args)
    return 1


if __name__ == "__main__":
    sys.exit(main())
