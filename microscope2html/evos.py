#!/usr/bin/env python3
"""Discovery and selection of EVOS M7000 scan directories.

An EVOS scan directory holds raw acquisition tiles named::

    scan_<slide>_R_p<pass>_<well>_A<area>f<field>d<channel>.TIF

plus, optionally, pre-stitched tile maps (``TM<N>``) exported by the
instrument and non-image artefacts (``.stitch`` files, scripts, assets).
This module classifies such files into *view units* — one viewer per
(slide, pass) — applying the project's channel policy.  Raw tiles are
always preferred; tile maps are only selected as a fallback for a unit
that has no raw tiles at all.
"""

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional

TILE_RE = re.compile(
    r"^scan_(?P<slide>.+)_(?P<type>R|TM\d*)_p(?P<pass>\d+)_"
    r"(?P<well>\d+)_A(?P<area>\d+)f(?P<field>\d+)d(?P<channel>\d+)"
    r"\.(?P<ext>[Tt][Ii][Ff][Ff]?)$"
)


@dataclass
class ViewUnit:
    """One viewer's worth of tiles: a slide acquired in one pass."""

    slide: str
    pas: str
    channel: Optional[int]
    raw_paths: List[str] = field(default_factory=list)
    tm_paths: List[str] = field(default_factory=list)
    tiles: List[tuple] = field(default_factory=list)

    @property
    def use_tm(self) -> bool:
        """True when the unit has no raw tiles and must use tile maps."""
        return not self.raw_paths

    @property
    def label(self) -> str:
        return f"slide '{self.slide}', pass {self.pas}"

    @property
    def slug(self) -> str:
        return slugify(self.slide)


def slugify(text: str) -> str:
    """Lowercase a slide name and turn runs of non-alphanumerics into '_'."""
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_") or "slide"


def parse_filename(name: str) -> Optional[Dict[str, str]]:
    """Parse an EVOS tile filename, or return None if it does not match."""
    match = TILE_RE.match(name)
    return match.groupdict() if match else None


def parse_channel(value) -> int:
    """Accept a channel as ``4`` or ``d4`` and return the integer."""
    text = str(value).strip()
    if text[:1] in ("d", "D"):
        text = text[1:]
    if not text.isdigit():
        raise ValueError(f"invalid channel {value!r} (expected e.g. 4 or d4)")
    return int(text)


def _fmt_channels(raw_by_channel: Dict[int, List[str]]) -> str:
    return ", ".join(f"d{channel}" for channel in sorted(raw_by_channel))


def discover_scan(directory: str, slide: Optional[str] = None,
                  channel: Optional[int] = None):
    """Classify an EVOS scan directory into view units.

    Args:
        directory: scan directory to inspect (top level only).
        slide: if given, keep only units from this exact slide name.
        channel: if given, the raw channel to use for every unit.

    Returns:
        ``(units, problems)`` — units ready for stitching and a list of
        human-readable problems for units that were skipped.
    """
    entries = []
    with os.scandir(directory) as it:
        for entry in it:
            if not entry.is_file():
                continue
            parsed = parse_filename(entry.name)
            if parsed:
                entries.append((entry.path, parsed))

    if not entries:
        return [], [
            f"no EVOS tiles found in {directory} (expected "
            f"scan_<slide>_R|TM<N>_p<pass>_<well>_A<area>f<field>d<channel>.TIF)"
        ]

    groups = {}
    for path, parsed in entries:
        key = (parsed["slide"], parsed["pass"])
        group = groups.setdefault(key, {"raw": {}, "tm": []})
        if parsed["type"] == "R":
            group["raw"].setdefault(int(parsed["channel"]), []).append(path)
        else:
            group["tm"].append(path)

    problems = []
    if slide is not None:
        matching = {key: value for key, value in groups.items() if key[0] == slide}
        if not matching:
            available = ", ".join(sorted({key[0] for key in groups})) or "none"
            return [], [
                f"no slide named {slide!r} in {directory} (available: {available})"
            ]
        groups = matching

    units = []
    for key in sorted(groups):
        slide_name, pas = key
        group = groups[key]
        raw_by_channel = group["raw"]
        tm_paths = sorted(group["tm"])
        label = f"slide '{slide_name}', pass {pas}"

        if not raw_by_channel:
            units.append(ViewUnit(slide_name, pas, None, [], tm_paths))
            continue

        if channel is not None:
            if channel not in raw_by_channel:
                problems.append(
                    f"{label}: requested channel d{channel} not found "
                    f"(available: {_fmt_channels(raw_by_channel)})"
                )
                continue
            selected = channel
        elif len(raw_by_channel) == 1:
            selected = next(iter(raw_by_channel))
        elif 4 in raw_by_channel:
            selected = 4
        else:
            problems.append(
                f"{label}: multiple channels ({_fmt_channels(raw_by_channel)}) "
                f"and no d4; re-run with --channel"
            )
            continue

        units.append(ViewUnit(slide_name, pas, selected,
                              sorted(raw_by_channel[selected]), tm_paths))

    return units, problems