"""
build_summary_tables.py

Turns results_summary_long.csv (produced by run_all_experiments.py) into the
paper-ready pivot tables:

  - table_twoarm.csv
      rows: instance id (1..5) + AVERAGE
      cols: H_{N}_obs_{p} for N in {20,30,40,50}, p in {3,4,5}   (12 columns)

  - table_multiarm_groups3.csv and table_multiarm_groups4.csv
      rows: instance id (1..5) + average
      cols: H_{N}_obs_{p} for N in {120,240,360}, p in {3,5}     (6 columns)

Each cell is the final best_fitness value from a single GA run (seed=123,
population=100, generations=200, crossover=1.0, mutation=0.14, tournament=2,
elitism=3), so the tables are fully reproducible from the supplement code.
"""

import argparse
from pathlib import Path

import pandas as pd

TWO_ARM_SUBJECTS = [20, 30, 40, 50]
TWO_ARM_P = [3, 4, 5]
MULTI_ARM_SUBJECTS = [120, 240, 360]
MULTI_ARM_P = [3, 5]
INSTANCE_IDS = [1, 2, 3, 4, 5]


def column_label(n: int, p: int) -> str:
    return f"H_{n}_obs_{p}"


def build_pivot(df: pd.DataFrame, subjects, ps, average_label: str) -> pd.DataFrame:
    columns = [column_label(n, p) for n in subjects for p in ps]
    table = pd.DataFrame(index=INSTANCE_IDS, columns=columns, dtype=float)

    for n in subjects:
        for p in ps:
            col = column_label(n, p)
            sub = df[(df["num_subjects"] == n) & (df["num_p"] == p)]
            for _, row in sub.iterrows():
                table.loc[row["instance_id"], col] = row["best_fitness"]

    table.index.name = "k"
    table.loc[average_label] = table.loc[INSTANCE_IDS].mean(numeric_only=True)
    return table


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--summary-csv", type=Path, default=Path("results/results_summary_long.csv"))
    ap.add_argument("--output-dir", type=Path, default=Path("results/tables"))
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.summary_csv)

    # --- Two-arm table ---
    two_arm_df = df[df["arm_type"] == "twoarm"]
    table_twoarm = build_pivot(two_arm_df, TWO_ARM_SUBJECTS, TWO_ARM_P, "AVERAGE")
    out = args.output_dir / "table_twoarm.csv"
    table_twoarm.to_csv(out)
    print(f"Wrote {out}")

    # --- Multi-arm tables (one per number of groups) ---
    multi_arm_df = df[df["arm_type"] == "multiarm"]
    for num_groups in sorted(multi_arm_df["num_groups"].unique()):
        sub = multi_arm_df[multi_arm_df["num_groups"] == num_groups]
        table = build_pivot(sub, MULTI_ARM_SUBJECTS, MULTI_ARM_P, "average")
        out = args.output_dir / f"table_multiarm_groups{int(num_groups)}.csv"
        table.to_csv(out)
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
