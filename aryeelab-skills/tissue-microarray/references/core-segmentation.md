# Core Segmentation

## Required Inputs

Minimum columns:
- `cell_id`
- `x_centroid`
- `y_centroid`
- `basic_qc`

Helpful optional columns:
- `transcript_counts`
- `ignore_from_core_analysis`
- existing analysis annotations that should be carried through

## Algorithm

1. Start from all cells not already marked ignored.
2. Infer column centers from x centroids using `kmeans` with the expected number of columns.
3. Assign each cell to the nearest inferred column center.
4. Infer row centers:
   - `global`: infer all row centers together from y centroids
   - `column-segments`: infer row centers separately within each column when columns are vertically offset
5. Build provisional `R#C#` labels from the inferred row and column identities.
6. For each provisional core:
   - use `basic_qc` cells only for the structural step
   - rasterize to `25 um` bins
   - binary close with `3x3`
   - label connected components
   - keep the largest connected component
   - compute its convex hull
   - keep every cell from that provisional core whose centroid lies inside the hull
7. Set cells outside the hull to missing final `core_id`.

## Interpretation

- `provisional_core_id` is the coarse grid assignment.
- `core_id` is the final convex-hull analysis assignment.
- `excluded_by_convex_hull_core == TRUE` means the cell was initially assigned to that core by the coarse grid but removed by the hull refinement.

## Defaults

- Default hull bin size: `25 um`
- Default morphology cleanup: `3x3` binary closing
- Default row mode: `global`

## Use Cases For Manual Overrides

- Ignore damaged strips with `--ignore-x-band`.
- Switch to `column-segments` when right and left blocks are vertically offset.
- Review cores with the lowest `hull_fraction` first.

## What Not To Do

- Do not use H&E to define the segmentation.
- Do not automatically reassign excluded cells to neighboring cores unless there is a separate validated reassignment procedure.
- Do not treat the convex hull as the true tissue edge; it is an analysis boundary.
