# GA optimal allocation — supplemental code

Code to reproduce the genetic algorithm allocation results reported in the
paper.

## Files

| File | Purpose |
|---|---|
| `ga_core.py` | GA logic: fitness, selection, crossover, mutation, I/O helpers. Supports any number of groups/arms. |
| `ga_cli.py` | Command-line entry point to reproduce **one** run on **one** instance. |
| `run_all_experiments.py` | Batch driver: runs the GA **once per instance**, over every two-arm and multi-arm instance used in the paper, with the fixed tuning parameters below. Writes `results_summary_long.csv`. |
| `build_summary_tables.py` | Turns `results_summary_long.csv` into the paper-ready pivot tables (`table_twoarm.csv`, `table_multiarm_groups3.csv`, `table_multiarm_groups4.csv`). |
| application | Folder with the code to evaluate the three-arm clinical trials for the weight maintenance trial. |

## Fixed tuning parameters (used for every single run)

```
population_size = 100
generations      = 200
elitism          = 3
tournament_size  = 2
mutation_rate    = 0.14
crossover_rate   = 1.00
seed             = 123
```

Because the seed is fixed and identical for every instance, each instance
is solved with a **single, reproducible GA run** — the same run used to
report the paper's results.

## Expected input layout

```
your_folder/
  covariate_matrices_twoarm/
    H_20_obs_3_id_1.csv ... H_50_obs_5_id_5.csv
  covariate_matrices_multiarm/
    H_120_obs_3_id_1.csv ... H_360_obs_5_id_5.csv
  moments_matrices/
    M_l_2_0.csv, M_l_2_1.csv, M_l_2_2.csv, M_l_4_0.csv
```

- **Two-arm** (`num_groups = 2`): N ∈ {20, 30, 40, 50}, p ∈ {3, 4, 5},
  id ∈ {1..5}. M file depends on p: `p=3 → M_l_2_0.csv`,
  `p=4 → M_l_2_1.csv`, `p=5 → M_l_2_2.csv`.
- **Multi-arm**: N ∈ {120, 240, 360}, p ∈ {3, 5}, id ∈ {1..5}. Every
  instance is solved **twice**: once with `num_groups = 3` and once with
  `num_groups = 4` (two independent sweeps over the full grid). M file
  depends on p: `p=3 → M_l_2_0.csv`, `p=5 → M_l_4_0.csv`.

By convention `H_*.csv` files have a leading subject-ID column (dropped
automatically) and `M_*.csv` files do not. If any of your files differ,
`ga_cli.py` exposes `--h-no-id-column` / `--m-has-id-column` flags, and
`ga_core.load_covariate_matrix` / `load_moment_matrix` accept the same
`id_column` argument directly.

## Usage

> **Note:** the commands below assume everything is downloaded into a
> **single folder**: the 4 `.py` files, `covariate_matrices_twoarm/`,
> `covariate_matrices_multiarm/`, and `moments_matrices/`, all at the same
> level:
>
> ```
> your_folder/
>     ga_core.py
>     ga_cli.py
>     run_all_experiments.py
>     build_summary_tables.py
>     covariate_matrices_twoarm/
>     covariate_matrices_multiarm/
>     moments_matrices/
> ```
>
> With this layout there is no need to type full paths: open a terminal,
> `cd` into that folder, and copy/paste the commands below as they are.
> (If your data is organized differently, replace the folder names in the
> commands with the corresponding full paths.)
>
> **Note on the `python` command:** depending on how Python is installed on
> your machine, the command to run it might be `python`, `python3`, `py`,
> or something else. Check which one works on your system by running
> `python --version` (and, if that fails, try `python3 --version` or
> `py --version`) in your terminal — whichever one returns a version number
> (instead of a "command not found" error) is the one to use. Replace
> `python` with that command in every command below.

1. Run every instance once (reproducible, seed=123 throughout):

   ```bash
   python run_all_experiments.py \
       --twoarm-dir covariate_matrices_twoarm \
       --multiarm-dir covariate_matrices_multiarm \
       --moments-dir moments_matrices \
       --output-root results
   ```

   This writes, per instance, a generation-history CSV and a subject
   assignment CSV under `results/twoarm/` or
   `results/multiarm/groups_{3,4}/`, plus one consolidated
   `results/results_summary_long.csv` (one row per instance: N, p, id,
   num_groups, best_fitness, elapsed time, file paths).

2. Build the paper tables:

   ```bash
   python build_summary_tables.py \
       --summary-csv results/results_summary_long.csv \
       --output-dir results/tables
   ```

   Produces:
   - `table_twoarm.csv` — rows `k = 1..5` + `AVERAGE`; columns
     `H_20_obs_3 ... H_50_obs_5` (12 columns, ordered by N then p).
   - `table_multiarm_groups3.csv` and `table_multiarm_groups4.csv` — rows
     `k = 1..5` + `average`; columns `H_120_obs_3 ... H_360_obs_5`
     (6 columns, ordered by N then p).

   Each cell is the final `best_fitness` from the single reproducible GA run for that
   instance.

3. (Optional) Reproduce a single instance manually, e.g. for a spot-check
   against the batch output:

   ```bash
   python ga_cli.py \
       --h-file covariate_matrices_twoarm/H_20_obs_3_id_1.csv \
       --m-file moments_matrices/M_l_2_0.csv \
       --num-subjects 20 --num-groups 2 --num-p 3 \
       --output-dir results/twoarm
   ```

4. (Application) Run the Weight Maintenance Diet Trial instance with **linear terms**,
   3 groups, p = 5, using
   `covariate_matrices_application/H_Weight_Maintenance_Diet_Trial.csv`
   and `moments_matrices/M_l_3_1.csv` (115 subjects):

   ```bash
   python ga_cli.py \
       --h-file covariate_matrices_application/H_Weight_Maintenance_Diet_Trial.csv \
       --m-file moments_matrices/M_l_3_1.csv \
       --num-subjects 115 --num-groups 3 --num-p 5 \
       --output-dir results/application
   ```

5. (Application) Run the expanded Weight Maintenance Diet Trial instance with **linear and quadratic terms**,
   3 groups, p = 8, using
   `covariate_matrices_application/H_Weight_Maintenance_Diet_Trial_expanded.csv`
   and `moments_matrices/M_q_4_1.csv` (115 subjects):

   ```bash
   python ga_cli.py \
       --h-file covariate_matrices_application/H_Weight_Maintenance_Diet_Trial_expanded.csv \
       --m-file moments_matrices/M_q_4_1.csv \
       --num-subjects 115 --num-groups 3 --num-p 8 \
       --output-dir results/application
   ```
