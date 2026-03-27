#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from scipy import ndimage as ndi
from scipy.spatial import ConvexHull
from sklearn.cluster import KMeans
from skimage import measure


def read_table(path):
    return pd.read_csv(path, compression="infer")


def write_table(df, path):
    path = Path(path)
    if path.suffix == ".gz":
        df.to_csv(path, index=False, compression="gzip")
    else:
        df.to_csv(path, index=False)


def infer_axis_clusters(values, n_centers, seed=1):
    km = KMeans(n_clusters=n_centers, random_state=seed, n_init=20)
    km.fit(np.asarray(values).reshape(-1, 1))
    centers = km.cluster_centers_.ravel()
    order = np.argsort(centers)
    remap = {int(order[i]): i + 1 for i in range(len(order))}
    labels = np.array([remap[int(x)] for x in km.labels_], dtype=int)
    return labels, np.sort(centers)


def assign_to_nearest_center(values, centers):
    values = np.asarray(values, dtype=float)
    centers = np.asarray(centers, dtype=float)
    return np.argmin(np.abs(values[:, None] - centers[None, :]), axis=1) + 1


def infer_column_y_segments(values, n_rows, bin_size, min_bin_cells, max_gap_bins):
    yb = np.floor(np.asarray(values, dtype=float) / bin_size).astype(int)
    bins = pd.Series(yb).value_counts().reset_index()
    bins.columns = ["yb", "n"]
    bins = bins.loc[bins["n"] >= min_bin_cells].sort_values("yb").reset_index(drop=True)
    if bins.empty:
        return np.quantile(values, np.linspace(0.1, 0.9, n_rows))
    gaps = bins["yb"].diff().fillna(0).astype(int)
    seg = (gaps > max_gap_bins).cumsum()
    segs = bins.groupby(seg, as_index=False)["yb"].mean().sort_values("yb")
    centers = segs["yb"].to_numpy() * bin_size + bin_size / 2
    if centers.shape[0] != n_rows:
        return np.quantile(values, np.linspace(0.1, 0.9, n_rows))
    return np.asarray(centers, dtype=float)


def convex_hull_polygon(points):
    if len(points) < 3:
        return None
    hull = ConvexHull(points)
    poly = points[hull.vertices]
    return np.vstack([poly, poly[0]])


def convex_hull_area(points):
    if len(points) < 3:
        return np.nan
    return float(ConvexHull(points).volume)


def connected_component_hull(frame, bin_size_um):
    basic = frame.loc[frame["basic_qc"].astype(bool)].copy()
    if basic.shape[0] < 3:
        points = frame[["x_centroid", "y_centroid"]].to_numpy()
        return np.ones(frame.shape[0], dtype=bool), {"n_components": 1, "largest_component_cells": int(basic.shape[0]), "second_component_cells": 0}, points

    x0 = float(basic["x_centroid"].min())
    y0 = float(basic["y_centroid"].min())
    xbin = np.floor((basic["x_centroid"].to_numpy() - x0) / bin_size_um).astype(int)
    ybin = np.floor((basic["y_centroid"].to_numpy() - y0) / bin_size_um).astype(int)
    grid = np.zeros((int(ybin.max()) + 2, int(xbin.max()) + 2), dtype=np.uint8)
    grid[ybin, xbin] = 1
    closed = ndi.binary_closing(grid, structure=np.ones((3, 3), dtype=bool))
    labels = measure.label(closed, connectivity=2)
    basic = basic.copy()
    basic["cc_label"] = labels[ybin, xbin]
    component_sizes = basic["cc_label"].value_counts().sort_values(ascending=False)
    keep_label = int(component_sizes.index[0])
    kept_basic = basic.loc[basic["cc_label"] == keep_label, ["x_centroid", "y_centroid"]].to_numpy()
    hull = convex_hull_polygon(kept_basic)
    if hull is None:
        keep_mask = np.ones(frame.shape[0], dtype=bool)
    else:
        path = MplPath(hull)
        points = frame[["x_centroid", "y_centroid"]].to_numpy()
        keep_mask = path.contains_points(points, radius=1.0)
    meta = {
        "n_components": int(component_sizes.shape[0]),
        "largest_component_cells": int(component_sizes.iloc[0]),
        "second_component_cells": int(component_sizes.iloc[1]) if component_sizes.shape[0] > 1 else 0,
    }
    return keep_mask, meta, kept_basic


def plot_spatial_core_map(df, out_path):
    sampled = df.sample(n=min(150000, df.shape[0]), random_state=1).copy()
    fig, ax = plt.subplots(figsize=(7.5, 10.5))
    keep = (~sampled["ignore_from_core_analysis"].astype(bool)) & sampled["core_id"].notna()
    if keep.any():
        ax.scatter(sampled.loc[keep, "x_centroid"], sampled.loc[keep, "y_centroid"], c=pd.factorize(sampled.loc[keep, "core_id"])[0], s=0.2, alpha=0.25, rasterized=True)
    drop = sampled["excluded_by_convex_hull_core"].astype(bool)
    if drop.any():
        ax.scatter(sampled.loc[drop, "x_centroid"], sampled.loc[drop, "y_centroid"], color="red", s=0.25, alpha=0.30, rasterized=True)
    ignored = sampled["ignore_from_core_analysis"].astype(bool)
    if ignored.any():
        ax.scatter(sampled.loc[ignored, "x_centroid"], sampled.loc[ignored, "y_centroid"], color="orange", s=0.25, alpha=0.30, rasterized=True)
    ax.set_title("Final TMA core calls")
    ax.set_xlabel("x centroid")
    ax.set_ylabel("y centroid")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_density_with_hulls(df, boundaries, out_path, density_bin_size, transcript_col):
    x = df["x_centroid"].to_numpy()
    y = df["y_centroid"].to_numpy()
    w = df[transcript_col].to_numpy()
    x_edges = np.arange(np.floor(x.min() / density_bin_size) * density_bin_size, np.ceil(x.max() / density_bin_size) * density_bin_size + density_bin_size, density_bin_size)
    y_edges = np.arange(np.floor(y.min() / density_bin_size) * density_bin_size, np.ceil(y.max() / density_bin_size) * density_bin_size + density_bin_size, density_bin_size)
    density, _, _ = np.histogram2d(y, x, bins=[y_edges, x_edges], weights=w)

    fig, ax = plt.subplots(figsize=(8.4, 10.8))
    im = ax.imshow(
        np.log1p(density),
        cmap="magma",
        origin="lower",
        extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
        interpolation="nearest",
    )
    for core_id, sub in boundaries.groupby("core_id"):
        sub = sub.sort_values("vertex_order")
        ax.plot(sub["x_um"], sub["y_um"], color="#7FDBFF", linewidth=1.2)
        ax.text(sub["x_um"].median(), sub["y_um"].median(), core_id, color="white", fontsize=6, ha="center", va="center")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(f"log1p summed {transcript_col} per {density_bin_size:.0f} um bin")
    ax.set_title("Convex-hull core boundaries on transcript density")
    ax.set_xlabel("x centroid")
    ax.set_ylabel("y centroid")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Infer TMA cores from scratch and refine them with convex hulls.")
    parser.add_argument("--cells", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--rows", type=int, required=True)
    parser.add_argument("--cols", type=int, required=True)
    parser.add_argument("--x-col", default="x_centroid")
    parser.add_argument("--y-col", default="y_centroid")
    parser.add_argument("--basic-qc-col", default="basic_qc")
    parser.add_argument("--transcript-col", default="transcript_counts")
    parser.add_argument("--ignore-col", default="ignore_from_core_analysis")
    parser.add_argument("--row-mode", choices=["global", "column-segments"], default="global")
    parser.add_argument("--bin-size-um", type=float, default=25.0)
    parser.add_argument("--column-segment-bin-size", type=float, default=100.0)
    parser.add_argument("--column-segment-min-bin-cells", type=int, default=10)
    parser.add_argument("--column-segment-max-gap-bins", type=int, default=2)
    parser.add_argument("--ignore-x-band", action="append", nargs=2, metavar=("START", "END"), help="Manual x-range to ignore before core calling.")
    parser.add_argument("--plot-density", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = read_table(args.cells).copy()
    required = ["cell_id", args.x_col, args.y_col, args.basic_qc_col]
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    rename_map = {
        args.x_col: "x_centroid",
        args.y_col: "y_centroid",
        args.basic_qc_col: "basic_qc",
    }
    if args.transcript_col in df.columns:
        rename_map[args.transcript_col] = "transcript_counts"
    if args.ignore_col in df.columns:
        rename_map[args.ignore_col] = "ignore_from_core_analysis"
    df = df.rename(columns=rename_map)
    if "ignore_from_core_analysis" not in df.columns:
        df["ignore_from_core_analysis"] = False

    if args.ignore_x_band:
        for start, end in args.ignore_x_band:
            start = float(start)
            end = float(end)
            df["ignore_from_core_analysis"] = df["ignore_from_core_analysis"].astype(bool) | ((df["x_centroid"] >= start) & (df["x_centroid"] <= end))

    ref = df.loc[~df["ignore_from_core_analysis"].astype(bool)].copy()
    _, x_centers = infer_axis_clusters(ref["x_centroid"], args.cols)
    ref["provisional_core_col"] = assign_to_nearest_center(ref["x_centroid"], x_centers)

    df["provisional_core_col"] = pd.NA
    df.loc[ref.index, "provisional_core_col"] = ref["provisional_core_col"].to_numpy()

    if args.row_mode == "global":
        _, y_centers = infer_axis_clusters(ref["y_centroid"], args.rows)
        df["provisional_core_row"] = pd.NA
        df.loc[ref.index, "provisional_core_row"] = assign_to_nearest_center(ref["y_centroid"], y_centers)
    else:
        df["provisional_core_row"] = pd.NA
        for cc in range(1, args.cols + 1):
            col_ref = ref.loc[ref["provisional_core_col"] == cc].copy()
            if col_ref.empty:
                continue
            y_centers = infer_column_y_segments(
                col_ref["y_centroid"].to_numpy(),
                args.rows,
                args.column_segment_bin_size,
                args.column_segment_min_bin_cells,
                args.column_segment_max_gap_bins,
            )
            df.loc[col_ref.index, "provisional_core_row"] = assign_to_nearest_center(col_ref["y_centroid"], y_centers)

    df["provisional_core_id"] = pd.NA
    keep = (~df["ignore_from_core_analysis"].astype(bool)) & df["provisional_core_row"].notna() & df["provisional_core_col"].notna()
    df.loc[keep, "provisional_core_id"] = [f"R{int(r)}C{int(c)}" for r, c in zip(df.loc[keep, "provisional_core_row"], df.loc[keep, "provisional_core_col"])]

    df["core_row"] = pd.NA
    df["core_col"] = pd.NA
    df["core_id"] = pd.NA
    df["excluded_by_convex_hull_core"] = False

    boundary_rows = []
    summary_rows = []
    analysis_mask = (~df["ignore_from_core_analysis"].astype(bool)) & df["provisional_core_id"].notna()

    for provisional_core_id, frame in df.loc[analysis_mask].groupby("provisional_core_id"):
        keep_mask, meta, kept_basic = connected_component_hull(frame, args.bin_size_um)
        keep_cell_ids = set(frame.loc[keep_mask, "cell_id"])
        idx = df["cell_id"].isin(frame["cell_id"])
        keep_idx = df["cell_id"].isin(keep_cell_ids)
        df.loc[idx & keep_idx, "core_row"] = frame["provisional_core_row"].iloc[0]
        df.loc[idx & keep_idx, "core_col"] = frame["provisional_core_col"].iloc[0]
        df.loc[idx & keep_idx, "core_id"] = provisional_core_id
        df.loc[idx & ~keep_idx, "excluded_by_convex_hull_core"] = True

        kept_points = frame.loc[keep_mask, ["x_centroid", "y_centroid"]].to_numpy()
        hull = convex_hull_polygon(kept_points)
        if hull is not None:
            for order_idx, (x_um, y_um) in enumerate(hull, start=1):
                boundary_rows.append({"core_id": provisional_core_id, "vertex_order": order_idx, "x_um": float(x_um), "y_um": float(y_um)})

        summary_rows.append(
            {
                "core_id": provisional_core_id,
                "provisional_core_row": int(frame["provisional_core_row"].iloc[0]),
                "provisional_core_col": int(frame["provisional_core_col"].iloc[0]),
                "old_n_cells": int(frame.shape[0]),
                "hull_n_cells": int(np.sum(keep_mask)),
                "excluded_n_cells": int(np.sum(~keep_mask)),
                "hull_fraction": float(np.mean(keep_mask)),
                "n_connected_components_basic_qc": meta["n_components"],
                "largest_basic_qc_component_cells": meta["largest_component_cells"],
                "second_basic_qc_component_cells": meta["second_component_cells"],
                "old_convex_hull_area_um2": convex_hull_area(frame[["x_centroid", "y_centroid"]].to_numpy()),
                "analysis_convex_hull_area_um2": convex_hull_area(kept_points),
            }
        )

    write_table(df, out_dir / "cell_annotations_with_core_calls.csv.gz")
    boundaries = pd.DataFrame(boundary_rows)
    boundaries.to_csv(out_dir / "convex_hull_boundaries.csv", index=False)
    pd.DataFrame(summary_rows).sort_values("hull_fraction").to_csv(out_dir / "per_core_hull_summary.csv", index=False)
    plot_spatial_core_map(df, out_dir / "spatial_core_map.png")

    if args.plot_density and "transcript_counts" in df.columns and not boundaries.empty:
        plot_density_with_hulls(
            df.loc[df["core_id"].notna()].copy(),
            boundaries,
            out_dir / "convex_hulls_on_transcript_density.png",
            density_bin_size=max(args.bin_size_um * 2, 50.0),
            transcript_col="transcript_counts",
        )


if __name__ == "__main__":
    main()
