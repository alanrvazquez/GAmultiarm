"""
run_all_experiments.py

Batch-runs the GA exactly ONCE per problem instance (fixed seed, therefore
deterministic and reproducible), using the fixed hyperparameters used
throughout the paper:

    population_size = 100
    generations      = 200
    elitism          = 3
    tournament_size  = 2
    mutation_rate    = 0.14
    crossover_rate   = 1.00
    seed             = 123   (identical for every single run)

Instance grids
--------------
Two-arm (num_groups = 2):
    N (num_subjects) in {20, 30, 40, 50}
    p (num_p)        in {3, 4, 5}
    id               in {1, 2, 3, 4, 5}
    H file: covariate_matrices_twoarm/H_{N}_obs_{p}_id_{id}.csv
    M file: moments_matrices/M_l_2_{p-3}.csv
            (p=3 -> M_l_2_0.csv, p=4 -> M_l_2_1.csv, p=5 -> M_l_2_2.csv)

Multi-arm (num_groups in {3, 4}, run as two separate sweeps over ALL
instances -- i.e. every (N, p, id) combination is solved once with 3 groups
and once with 4 groups):
    N (num_subjects) in {120, 240, 360}
    p (num_p)        in {3, 5}
    id               in {1, 2, 3, 4, 5}
    H file: covariate_matrices_multiarm/H_{N}_obs_{p}_id_{id}.csv
    M file: moments_matrices/M_l_2_0.csv (p=3) or moments_matrices/M_l_4_0.csv (p=5)

All raw per-instance outputs (generation history + final assignment) are
written under --output-root, and a single long-format CSV
(results_summary_long.csv) with one row per instance is written as results
come in (so an interrupted run still leaves usable partial output).
build_summary_tables.py turns that file into the paper-ready pivot tables.
"""

import argparse
import csv
from pathlib import Path
from typing import Tuple

from ga_core import load_covariate_matrix, load_moment_matrix, run_single_instance

# Fixed GA hyperparameters used for every instance in the paper.
DEFAULT_PARAMS = dict(
    population_size=100,
    generations=200,
    crossover_rate=1.0,
    mutation_rate=0.14,
    tournament_size=2,
    elitism=3,
    seed=123,
)

TWO_ARM_SUBJECTS = [20, 30, 40, 50]
TWO_ARM_P = [3, 4, 5]
MULTI_ARM_SUBJECTS = [120, 240, 360]
MULTI_ARM_P = [3, 5]
INSTANCE_IDS = [1, 2, 3, 4, 5]
MULTI_ARM_GROUP_COUNTS = [3, 4]

TWO_ARM_M_FILE = {3: "M_l_2_0.csv", 4: "M_l_2_1.csv", 5: "M_l_2_2.csv"}
MULTI_ARM_M_FILE = {3: "M_l_2_0.csv", 5: "M_l_4_0.csv"}


def h_filename(n: int, p: int, k: int) -> str:
    return f"H_{n}_obs_{p}_id_{k}.csv"


def build_two_arm_instances(covariate_dir: Path, moments_dir: Path):
    instances = []
    for n in TWO_ARM_SUBJECTS:
        for p in TWO_ARM_P:
            for k in INSTANCE_IDS:
                instances.append(dict(
                    arm_type="twoarm",
                    num_groups=2,
                    num_subjects=n,
                    num_p=p,
                    instance_id=k,
                    h_path=covariate_dir / h_filename(n, p, k),
                    m_path=moments_dir / TWO_ARM_M_FILE[p],
                ))
    return instances


def build_multi_arm_instances(covariate_dir: Path, moments_dir: Path):
    instances = []
    for num_groups in MULTI_ARM_GROUP_COUNTS:
        for n in MULTI_ARM_SUBJECTS:
            for p in MULTI_ARM_P:
                for k in INSTANCE_IDS:
                    instances.append(dict(
                        arm_type="multiarm",
                        num_groups=num_groups,
                        num_subjects=n,
                        num_p=p,
                        instance_id=k,
                        h_path=covariate_dir / h_filename(n, p, k),
                        m_path=moments_dir / MULTI_ARM_M_FILE[p],
                    ))
    return instances


def run_all(twoarm_dir: Path, multiarm_dir: Path, moments_dir: Path, output_root: Path, params: dict):
    instances = (
        build_two_arm_instances(twoarm_dir, moments_dir)
        + build_multi_arm_instances(multiarm_dir, moments_dir)
    )

    output_root.mkdir(parents=True, exist_ok=True)
    summary_path = output_root / "results_summary_long.csv"
    summary_rows = []

    with open(summary_path, "w", newline="", encoding="utf-8") as f_out:
        writer = None

        for inst in instances:
            h_path, m_path = inst["h_path"], inst["m_path"]

            if not h_path.exists():
                print(f"[SKIP] missing H file: {h_path}")
                continue
            if not m_path.exists():
                print(f"[SKIP] missing M file: {m_path}")
                continue

            H, _ = load_covariate_matrix(h_path)
            M = load_moment_matrix(m_path)

            if inst["arm_type"] == "twoarm":
                out_dir = output_root / "twoarm"
            else:
                out_dir = output_root / "multiarm" / f"groups_{inst['num_groups']}"

            result = run_single_instance(
                H=H,
                M=M,
                instance_stem=h_path.stem,
                num_subjects=inst["num_subjects"],
                num_groups=inst["num_groups"],
                num_p=inst["num_p"],
                output_dir=out_dir,
                **params,
            )
            result["arm_type"] = inst["arm_type"]
            result["instance_id"] = inst["instance_id"]
            result["h_file"] = str(h_path)
            result["m_file"] = str(m_path)

            if writer is None:
                writer = csv.DictWriter(f_out, fieldnames=list(result.keys()))
                writer.writeheader()
            writer.writerow(result)
            f_out.flush()  # incremental write: partial runs stay usable

            summary_rows.append(result)

    print(f"\nWrote {len(summary_rows)} instance results to {summary_path}")
    return summary_path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--covariate-root", type=Path, default=None,
                    help="Convenience option: a single directory containing both "
                         "covariate_matrices_twoarm/ and covariate_matrices_multiarm/ "
                         "as subfolders. Ignored if --twoarm-dir / --multiarm-dir are given.")
    p.add_argument("--twoarm-dir", type=Path, default=None,
                    help="Directory holding the two-arm H_*.csv files directly "
                         "(e.g. .../covariate_matrices_twoarm). Overrides --covariate-root.")
    p.add_argument("--multiarm-dir", type=Path, default=None,
                    help="Directory holding the multi-arm H_*.csv files directly "
                         "(e.g. .../covariate_matrices_multiarm). Overrides --covariate-root.")
    p.add_argument("--moments-dir", type=Path, default=Path("moments_matrices"))
    p.add_argument("--output-root", type=Path, default=Path("results"))
    p.add_argument("--population-size", type=int, default=DEFAULT_PARAMS["population_size"])
    p.add_argument("--generations", type=int, default=DEFAULT_PARAMS["generations"])
    p.add_argument("--crossover-rate", type=float, default=DEFAULT_PARAMS["crossover_rate"])
    p.add_argument("--mutation-rate", type=float, default=DEFAULT_PARAMS["mutation_rate"])
    p.add_argument("--tournament-size", type=int, default=DEFAULT_PARAMS["tournament_size"])
    p.add_argument("--elitism", type=int, default=DEFAULT_PARAMS["elitism"])
    p.add_argument("--seed", type=int, default=DEFAULT_PARAMS["seed"])
    return p.parse_args()


def resolve_arm_dirs(args) -> Tuple[Path, Path]:
    """
    Resolve the two-arm / multi-arm covariate directories from either the
    explicit --twoarm-dir / --multiarm-dir flags, or the convenience
    --covariate-root flag (root/covariate_matrices_twoarm, root/covariate_matrices_multiarm).
    """
    if args.twoarm_dir is not None or args.multiarm_dir is not None:
        if args.twoarm_dir is None or args.multiarm_dir is None:
            raise SystemExit("Pass BOTH --twoarm-dir and --multiarm-dir together "
                              "(or use --covariate-root instead).")
        return args.twoarm_dir, args.multiarm_dir

    if args.covariate_root is not None:
        return (args.covariate_root / "covariate_matrices_twoarm",
                args.covariate_root / "covariate_matrices_multiarm")

    raise SystemExit("Provide either --covariate-root, or both --twoarm-dir and --multiarm-dir.")


def main():
    args = parse_args()
    twoarm_dir, multiarm_dir = resolve_arm_dirs(args)
    params = dict(
        population_size=args.population_size,
        generations=args.generations,
        crossover_rate=args.crossover_rate,
        mutation_rate=args.mutation_rate,
        tournament_size=args.tournament_size,
        elitism=args.elitism,
        seed=args.seed,
    )
    run_all(twoarm_dir, multiarm_dir, args.moments_dir, args.output_root, params)


if __name__ == "__main__":
    main()
