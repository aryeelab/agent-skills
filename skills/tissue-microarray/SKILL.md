---
name: tissue-microarray
description: Segment tissue microarray cores from raw cell centroid tables, then refine each inferred core with a largest-connected-component plus convex-hull boundary. Use when Codex needs to infer TMA core identities, clean up spill-in between neighboring cores, generate core boundary overlays, or propagate final core labels into downstream spatial transcriptomics tables. Best for Xenium-like or similar per-cell tables with x/y centroids, QC flags, and an expected TMA grid such as 7x4.
---

# Tissue Microarray

Use this skill to infer TMA cores directly from cell centroids and convert them into convex-hull analysis boundaries.

## Workflow

1. Start from a per-cell table with at least:
   - `cell_id`
   - `x_centroid`
   - `y_centroid`
   - `basic_qc`
2. Decide the expected TMA grid, usually `--rows` and `--cols`.
3. If needed, define manual ignore bands before core inference, for example a damaged middle strip.
4. Infer a coarse grid from centroids:
   - infer columns from x positions
   - infer rows globally or within each column
5. Within each inferred provisional core:
   - use `basic_qc` cells to define tissue structure
   - rasterize centroids to `25 um` bins unless there is a strong reason to change scale
   - apply a `3x3` binary closing
   - label connected components
   - keep the largest connected component
   - compute its convex hull
   - define the final analysis core as all cells from that provisional core whose centroids fall inside the hull
6. Use H&E only post-hoc as QC, never to define the core assignment.

## Main Script

Run the segmentation with:

```bash
python scripts/segment_tma.py \
  --cells /path/to/cells.csv.gz \
  --out-dir /path/to/output_dir \
  --rows 7 \
  --cols 4 \
  --plot-density
```

Optional manual ignore bands:

```bash
python scripts/segment_tma.py \
  --cells /path/to/cells.csv.gz \
  --out-dir /path/to/output_dir \
  --rows 7 \
  --cols 4 \
  --ignore-x-band 5200 6350
```

If rows are vertically offset across columns, use:

```bash
python scripts/segment_tma.py \
  --cells /path/to/cells.csv.gz \
  --out-dir /path/to/output_dir \
  --rows 7 \
  --cols 4 \
  --row-mode column-segments
```

This writes:
- `cell_annotations_with_core_calls.csv.gz`
- `per_core_hull_summary.csv`
- `convex_hull_boundaries.csv`
- `spatial_core_map.png`
- `convex_hulls_on_transcript_density.png` if `transcript_counts` exists and `--plot-density` is set

## Propagating Labels

If another table should inherit the final core labels by `cell_id`, use:

```bash
python scripts/propagate_core_labels.py \
  --source /path/to/output_dir/cell_annotations_with_core_calls.csv.gz \
  --target /path/to/other_table.csv.gz \
  --out /path/to/other_table_with_core_labels.csv.gz
```

## QC Rules

- Inspect `per_core_hull_summary.csv` first.
- Low `hull_fraction` is the main signal that the coarse grid assignment was pulling in spill-in.
- Convex hulls are acceptable as both the displayed and analysis boundary even if they contain some empty space.
- Use H&E afterward only to ask whether the hulls land on plausible tissue islands.

## Resources

- Read [references/core-segmentation.md](references/core-segmentation.md) for the exact algorithm and caveats.
- Use [scripts/segment_tma.py](scripts/segment_tma.py) for end-to-end segmentation.
- Use [scripts/propagate_core_labels.py](scripts/propagate_core_labels.py) to copy final core labels into downstream tables.
