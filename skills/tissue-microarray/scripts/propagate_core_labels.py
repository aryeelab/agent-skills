#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def read_table(path):
    return pd.read_csv(path, compression="infer")


def write_table(df, path):
    path = Path(path)
    if path.suffix == ".gz":
        df.to_csv(path, index=False, compression="gzip")
    else:
        df.to_csv(path, index=False)


def parse_args():
    parser = argparse.ArgumentParser(description="Copy final TMA core labels by cell_id into another table.")
    parser.add_argument("--source", required=True, help="cell_annotations_with_core_calls.csv.gz from segment_tma_from_scratch.py")
    parser.add_argument("--target", required=True, help="Table to receive the final core columns.")
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    source = read_table(args.source)
    target = read_table(args.target)

    keep_cols = [
        "cell_id",
        "provisional_core_row",
        "provisional_core_col",
        "provisional_core_id",
        "core_row",
        "core_col",
        "core_id",
        "excluded_by_convex_hull_core",
        "ignore_from_core_analysis",
    ]
    missing = [col for col in keep_cols if col not in source.columns]
    if missing:
        raise ValueError(f"Source table is missing required propagated columns: {missing}")

    propagated = source[keep_cols].drop_duplicates("cell_id")
    target = target.drop(columns=[col for col in keep_cols if col in target.columns and col != "cell_id"], errors="ignore")
    out = target.merge(propagated, on="cell_id", how="left")
    write_table(out, args.out)


if __name__ == "__main__":
    main()
