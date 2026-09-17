"""
ga_core.py

Genetic algorithm for optimal subject-to-group (arm) allocation under a
minimax average GSP variance criterion.

This module is an English-language, reorganized version of the original
allocation script. The optimization logic itself (selection, crossover,
mutation, fitness) is unchanged. It has been split into reusable
functions/classes so that many problem instances can be solved
programmatically (see run_all_experiments.py) instead of one at a time from
the command line, and so that the number of groups/arms is no longer
hard-coded to two.
"""

import csv
import random
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
@dataclass
class GAConfig:
    """Genetic algorithm configuration for one problem instance."""

    # Problem size
    num_subjects: int   # total number of subjects available
    num_groups: int      # number of arms/groups to allocate subjects into
    num_p: int            # number of covariates (columns of H / size of M)

    # Chromosome size 
    # If chromosome_length == num_subjects, every chromosome is a permutation
    # of the same subject set, i.e. every subject is allocated to a group.
    chromosome_length: int

    # GA hyperparameters
    population_size: int
    generations: int
    crossover_rate: float
    mutation_rate: float
    tournament_size: int
    elitism: int
    seed: int


@dataclass
class ProblemData:
    """Precomputed, instance-specific data shared by all fitness evaluations."""

    H: np.ndarray                              # (num_subjects, num_p) covariate matrix
    M: np.ndarray                               # (num_p, num_p) moment matrix
    group_sizes: List[int]                      # sizes of the groups; sum == chromosome_length
    hi_hit: List[np.ndarray] = field(default_factory=list)  # h_i h_i^T per subject

    def __post_init__(self):
        if not self.hi_hit:
            self.hi_hit = [
                self.H[i, :].reshape(-1, 1) @ self.H[i, :].reshape(1, -1)
                for i in range(self.H.shape[0])
            ]


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def balanced_group_sizes(num_subjects: int, num_groups: int) -> List[int]:
    """
    Split num_subjects as evenly as possible into num_groups groups.
    Any remainder is distributed to the last groups (largest-remainder rule),
    e.g. 19 subjects / 4 groups -> [4, 5, 5, 5].
    """
    base, remainder = divmod(num_subjects, num_groups)
    sizes = [base] * num_groups
    for i in range(num_groups - remainder, num_groups):
        sizes[i] += 1
    return sizes


def load_covariate_matrix(path: Path, id_column: bool = True) -> Tuple[np.ndarray, pd.DataFrame]:
    """
    Load an H (covariate) matrix from CSV.
    By convention these files store a subject ID in the first column,
    followed by the numeric covariates; pass id_column=False if a given
    file has no such leading column.
    """
    df = pd.read_csv(path)
    H = df.iloc[:, 1:].to_numpy() if id_column else df.to_numpy()
    return H.astype(float), df


def load_moment_matrix(path: Path, id_column: bool = False) -> np.ndarray:
    """
    Load an M (moment) matrix from CSV. By convention these files contain
    only the p x p numeric matrix (no ID column). Pass id_column=True if a
    given file has a leading ID/index column that must be dropped.
    """
    df = pd.read_csv(path)
    M = df.iloc[:, 1:].to_numpy() if id_column else df.to_numpy()
    return M.astype(float)


# --------------------------------------------------------------------------
# Fitness functions
# --------------------------------------------------------------------------
def fitness(chromosome: List[int], cfg: GAConfig, data: ProblemData) -> float:
    """
    Minimax average GSP criterion.

    The chromosome (a permutation of subject indices) is split, in order,
    into cfg.num_groups consecutive blocks whose sizes are data.group_sizes.
    For each group k, the information matrix is the sum of h_i h_i^T over the
    subjects assigned to that group; the criterion value for the group is
    trace(inverse(information matrix) @ M). The chromosome's fitness is the
    WORST (maximum) value across groups, which run_ga MINIMIZES (a minimax
    design criterion). Singular information matrices are penalized with a
    large constant instead of raising an error.
    """
    BIG_VALUE = 1e15
    group_weights = [1.0] * cfg.num_groups  # equal weighting across groups

    var_k = [0.0] * cfg.num_groups
    start = 0
    for k, size in enumerate(data.group_sizes):
        sum_xhht = np.zeros((cfg.num_p, cfg.num_p))
        for i in range(start, start + size):
            subject_i = chromosome[i]
            sum_xhht += data.hi_hit[subject_i]
        start += size

        try:
            inv_sum = np.linalg.inv(sum_xhht)
            var_k[k] = group_weights[k] * np.trace(inv_sum @ data.M)
        except np.linalg.LinAlgError:
            var_k[k] = BIG_VALUE

    return max(var_k)


def fitness_2arm_pairwise(chromosome: List[int], cfg: GAConfig, H: np.ndarray, M: np.ndarray) -> float:
    """
    Alternative, pairwise-sum formulation for the balanced 2-arm case.
    Kept for reference only: the main pipeline (run_all_experiments.py)
    uses `fitness` for every instance type, two-arm and multi-arm alike.
    Requires len(chromosome) to be divisible by cfg.num_groups.
    """
    n = len(chromosome)
    num_parts = cfg.num_groups
    if n % num_parts != 0:
        raise ValueError("chromosome length must be divisible by cfg.num_groups.")

    part_length = n // num_parts
    hth_inv = np.linalg.inv(H.T @ H)
    obj = 0.0

    for i in range(n - 1):
        idx_i = chromosome[i]
        part_of_i = i // part_length
        Ji = hth_inv @ H[idx_i, :, None] @ H[idx_i, :, None].T

        for j in range(i + 1, n):
            idx_j = chromosome[j]
            part_of_j = j // part_length
            Jj = hth_inv @ H[idx_j, :, None] @ H[idx_j, :, None].T
            Cij = np.trace(Ji @ Jj @ hth_inv @ M)
            obj += (2 * part_of_i - 1) * (2 * part_of_j - 1) * Cij

    return 1.0 / (1.0 + obj)


# --------------------------------------------------------------------------
# GA operators
# --------------------------------------------------------------------------
def create_chromosome(cfg: GAConfig) -> List[int]:
    return random.sample(range(cfg.num_subjects), cfg.chromosome_length)


def create_population(cfg: GAConfig) -> List[List[int]]:
    return [create_chromosome(cfg) for _ in range(cfg.population_size)]


def tournament_selection(population: List[List[int]], cfg: GAConfig, data: ProblemData) -> List[int]:
    contenders = random.sample(population, cfg.tournament_size)
    return min(contenders, key=lambda ch: fitness(ch, cfg, data))  # minimization


def order_crossover(p1: List[int], p2: List[int]) -> Tuple[List[int], List[int]]:
    """Order crossover (OX) for permutations; preserves gene uniqueness."""
    n = len(p1)
    a, b = sorted(random.sample(range(n), 2))

    def make_child(pa, pb):
        child = [-1] * n
        child[a:b + 1] = pa[a:b + 1]
        fill_values = [g for g in pb if g not in child]
        idx = 0
        for i in range(n):
            if child[i] == -1:
                child[i] = fill_values[idx]
                idx += 1
        return child

    return make_child(p1, p2), make_child(p2, p1)


def swap_mutation(ch: List[int]) -> List[int]:
    """Swap two positions; preserves uniqueness and the valid gene range."""
    c = ch[:]
    i, j = random.sample(range(len(c)), 2)
    c[i], c[j] = c[j], c[i]
    return c


def run_ga(cfg: GAConfig, data: ProblemData) -> Tuple[List[int], float, List[float]]:
    if cfg.chromosome_length > cfg.num_subjects:
        raise ValueError("chromosome_length must be <= num_subjects")

    random.seed(cfg.seed)
    population = create_population(cfg)
    history_best = []

    for gen in range(cfg.generations):
        population.sort(key=lambda ch: fitness(ch, cfg, data))  # ascending (minimize)
        next_pop = population[:cfg.elitism]

        while len(next_pop) < cfg.population_size:
            p1 = tournament_selection(population, cfg, data)
            p2 = tournament_selection(population, cfg, data)

            if random.random() < cfg.crossover_rate:
                c1, c2 = order_crossover(p1, p2)
            else:
                c1, c2 = p1[:], p2[:]

            if random.random() < cfg.mutation_rate:
                c1 = swap_mutation(c1)
            if random.random() < cfg.mutation_rate:
                c2 = swap_mutation(c2)

            next_pop.append(c1)
            if len(next_pop) < cfg.population_size:
                next_pop.append(c2)

        population = next_pop
        best_now = min(population, key=lambda ch: fitness(ch, cfg, data))
        best_fitness = fitness(best_now, cfg, data)
        history_best.append(best_fitness)

        if gen % 10 == 0:
            print(f"  Gen {gen}: best={best_fitness:.6f}")

    best = min(population, key=lambda ch: fitness(ch, cfg, data))
    return best, fitness(best, cfg, data), history_best


# --------------------------------------------------------------------------
# Output helpers
# --------------------------------------------------------------------------
def save_assignment(chromosome: List[int], group_sizes: List[int], out_file: Path) -> Path:
    """
    Save the subject -> group assignment implied by `chromosome` (split into
    consecutive blocks of size group_sizes[0], group_sizes[1], ...).

    Generalized to an arbitrary number of groups. For two groups this
    function still reports the -1/+1 code (group_code_pm1) for backward
    compatibility, in addition to the 1-based group_index that works for
    any number of arms.
    """
    num_subj = sum(group_sizes)
    if len(chromosome) < num_subj:
        raise ValueError("chromosome is shorter than sum(group_sizes)")

    group_index: List[Optional[int]] = [None] * len(chromosome)
    start = 0
    for g, size in enumerate(group_sizes):
        for i in range(start, start + size):
            subj = chromosome[i]
            group_index[subj] = g + 1  # 1-based group label
        start += size

    two_arm = len(group_sizes) == 2

    with open(out_file, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        header = ["subject_id", "group_index"] + (["group_code_pm1"] if two_arm else [])
        w.writerow(header)
        for subj_id, g in enumerate(group_index):
            if g is None:
                continue  # subject not used (chromosome_length < num_subjects)
            row = [subj_id, g]
            if two_arm:
                row.append(-1 if g == 1 else 1)
            w.writerow(row)

    return out_file


def save_generation_history(history: List[float], elapsed: float, group_sizes: List[int],
                             best_chromosome: List[int], out_file: Path) -> Path:
    """Save per-generation best fitness plus the group split of the final best chromosome."""
    df = pd.DataFrame({"iter": range(1, len(history) + 1), "best": history})
    df["elapsed"] = pd.NA
    df.loc[0, "elapsed"] = elapsed

    chrom = np.asarray(best_chromosome).ravel()
    starts = [0]
    for s in group_sizes[:-1]:
        starts.append(starts[-1] + s)

    groups = {
        f"group_{g + 1}": chrom[starts[g]: starts[g] + group_sizes[g]]
        for g in range(len(group_sizes))
    }
    groups_df = pd.DataFrame({k: pd.Series(v) for k, v in groups.items()})

    max_rows = max(len(df), len(groups_df))
    out_df = pd.concat(
        [df.reindex(range(max_rows)), groups_df.reindex(range(max_rows))], axis=1
    )
    out_df.to_csv(out_file, index=False)
    return out_file


# --------------------------------------------------------------------------
# One-instance orchestration
# --------------------------------------------------------------------------
def run_single_instance(
    H: np.ndarray,
    M: np.ndarray,
    instance_stem: str,
    num_subjects: int,
    num_groups: int,
    num_p: int,
    population_size: int,
    generations: int,
    crossover_rate: float,
    mutation_rate: float,
    tournament_size: int,
    elitism: int,
    seed: int,
    output_dir: Path,
    verbose: bool = True,
) -> dict:
    """
    Run the GA once (deterministic, seeded) on one problem instance and write
    the per-generation history and the final subject assignment to
    output_dir. Returns a summary dict; run_all_experiments.py collects one
    of these per instance to build the paper tables.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cfg = GAConfig(
        num_subjects=num_subjects,
        num_groups=num_groups,
        num_p=num_p,
        chromosome_length=num_subjects,
        population_size=population_size,
        generations=generations,
        crossover_rate=crossover_rate,
        mutation_rate=mutation_rate,
        tournament_size=tournament_size,
        elitism=elitism,
        seed=seed,
    )
    group_sizes = balanced_group_sizes(num_subjects, num_groups)
    data = ProblemData(H=H, M=M, group_sizes=group_sizes)

    if verbose:
        print(f"[{instance_stem}] N={num_subjects}, groups={num_groups}, p={num_p}, "
              f"group_sizes={group_sizes}")

    start = time.perf_counter()
    best_chromosome, best_fit, history = run_ga(cfg, data)
    elapsed = time.perf_counter() - start

    c = f"{cfg.crossover_rate:.2f}".replace(".", "p")
    m = f"{cfg.mutation_rate:.2f}".replace(".", "p")
    tag = f"P{cfg.population_size}_G{cfg.generations}_C{c}_M{m}_T{cfg.tournament_size}_E{cfg.elitism}"

    results_file = output_dir / f"GA_results_{instance_stem}_{tag}.csv"
    assignment_file = output_dir / f"GA_assignment_{instance_stem}_{tag}.csv"

    save_generation_history(history, elapsed, group_sizes, best_chromosome, results_file)
    save_assignment(best_chromosome, group_sizes, assignment_file)

    if verbose:
        print(f"[{instance_stem}] best_fitness={best_fit:.6f}  elapsed={elapsed:.2f}s")

    return {
        "instance_stem": instance_stem,
        "num_subjects": num_subjects,
        "num_groups": num_groups,
        "num_p": num_p,
        "population_size": population_size,
        "generations": generations,
        "crossover_rate": crossover_rate,
        "mutation_rate": mutation_rate,
        "tournament_size": tournament_size,
        "elitism": elitism,
        "seed": seed,
        "best_fitness": best_fit,
        "elapsed_seconds": elapsed,
        "results_file": str(results_file),
        "assignment_file": str(assignment_file),
    }
