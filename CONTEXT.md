# CONTEXT.md — microscope2html glossary

## Slide
A physical glass slide with tissue on it, placed under the microscope. All input TIFF files
for a single invocation belong to the same slide.

## Tile
A single TIFF file captured by the microscope. Contains one image frame with OME-TIFF
metadata including stage coordinates (`PositionX`, `PositionY`) and pixel size (`PhysicalSizeX`).

## ROI (Region of Interest)
A cluster of spatially adjacent tiles from the same slide. ROIs are detected automatically
by distance-based clustering of tile stage coordinates. Tiles within the same ROI are close
enough to be part of a single contiguous scan; tiles from different ROIs are separated by
substantial empty space.

An ROI may be:
- **Multi-tile**: a grid of overlapping tiles that form a contiguous image (stitched by Fiji).
- **Single-tile**: a lone tile with no neighbours (passed through without stitching).

## Cluster threshold
The configurable multiplier (N× tile diagonal) used by distance-based clustering to group
tiles into ROIs. Tiles whose centre-to-centre distance is less than this threshold are
assigned to the same ROI.

## Margin
The configurable pixel gap inserted between ROI bounding boxes when collapsing empty bands
on the compressed canvas. Prevents ROIs from visually touching.

## Empty band
A contiguous horizontal or vertical span of the coordinate space that contains no tile data
from any ROI. Empty bands are collapsed (removed) when laying out ROIs in the viewer.

## Compressed canvas
The coordinate space of the OpenSeadragon viewer after empty bands have been collapsed
and configurable margins have been inserted between ROIs. ROIs retain their relative
directional ordering (e.g., ROI-A that was above-and-left of ROI-B remains above-and-left),
but the empty distance between them is eliminated.

## Stitched image
The output of Fiji's Grid/Collection stitching plugin for a single ROI — a single TIFF
where tiles have been placed at their relative coordinates and overlapping edges blended.

## DZI (Deep Zoom Image)
A multi-resolution tile pyramid generated from a stitched image. Each ROI gets its own
DZI, embedded as base64-encoded JPEG tiles in the output HTML. OpenSeadragon loads the
appropriate resolution tiles based on the current zoom level.