"""
ga_cli.py

Command-line entry point to reproduce ONE run of the genetic algorithm on
ONE problem instance (one H file + one M file). This mirrors the original
standalone script, but with named arguments instead of positional ones for
clarity in the supplement.

Example (two-arm instance):
    python ga_cli.py \
        --h-file covariate_matrices_twoarm/H_20_obs_3_id_1.csv \
        --m-file moments_matrices/M_l_2_0.csv \
        --num-subjects 20 --num-groups 2 --num-p 3 \
        --population-size 100 --generations 200 \
        --crossover-rate 1.0 --mutation-rate 0.14 \
        --tournament-size 2 --elitism 3 --seed 123 \
        --output-dir results/twoarm

Example (multi-arm instance, 3 groups):
    python ga_cli.py \
        --h-file covariate_matrices_multiarm/H_120_obs_3_id_1.csv \
        --m-file moments_matrices/M_l_2_0.csv \
        --num-subjects 120 --num-groups 3 --num-p 3 \
        --population-size 100 --generations 200 \
        --crossover-rate 1.0 --mutation-rate 0.14 \
        --tournament-size 2 --elitism 3 --seed 123 \
        --output-dir results/multiarm/groups_3
"""

import argparse
from pathlib import Path

from ga_core import load_covariate_matrix, load_moment_matrix, run_single_instance


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--h-file", required=True, type=Path)
    p.add_argument("--m-file", required=True, type=Path)
    p.add_argument("--num-subjects", required=True, type=int)
    p.add_argument("--num-groups", required=True, type=int)
    p.add_argument("--num-p", required=True, type=int)
    p.add_argument("--population-size", type=int, default=100)
    p.add_argument("--generations", type=int, default=200)
    p.add_argument("--crossover-rate", type=float, default=1.0)
    p.add_argument("--mutation-rate", type=float, default=0.14)
    p.add_argument("--tournament-size", type=int, default=2)
    p.add_argument("--elitism", type=int, default=3)
    p.add_argument("--seed", type=int, default=123)
    p.add_argument("--output-dir", type=Path, default=Path("."))
    p.add_argument("--h-no-id-column", action="store_true",
                    help="Pass this if the H file has NO leading ID column.")
    p.add_argument("--m-has-id-column", action="store_true",
                    help="Pass this if the M file HAS a leading ID/index column.")
    return p.parse_args()


def main():
    args = parse_args()

    H, _ = load_covariate_matrix(args.h_file, id_column=not args.h_no_id_column)
    M = load_moment_matrix(args.m_file, id_column=args.m_has_id_column)

    summary = run_single_instance(
        H=H,
        M=M,
        instance_stem=args.h_file.stem,
        num_subjects=args.num_subjects,
        num_groups=args.num_groups,
        num_p=args.num_p,
        population_size=args.population_size,
        generations=args.generations,
        crossover_rate=args.crossover_rate,
        mutation_rate=args.mutation_rate,
        tournament_size=args.tournament_size,
        elitism=args.elitism,
        seed=args.seed,
        output_dir=args.output_dir,
    )

    print("\n=== Summary ===")
    for k, v in summary.items():
        print(f"{k}: {v}")


if __name__ == "__main__":
    main()
