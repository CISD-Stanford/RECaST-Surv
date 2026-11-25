# RECaST-Surv

This repository contains the RECaST-Surv implementation and simulation code accompanying “RECaST-Surv: A Calibrated Borrowing Method for Survival Endpoints in Unequal Randomized Trials.” It includes the main simulation settings, robustness scenarios H1–H5, and sensitivity analyses A–D, together with scripts for generating data, running analyses, and summarizing newly generated results. See the manuscript and appendices for the statistical definitions, assumptions, and interpretation.

This guide follows the manuscript version dated November 24, 2026. Run the commands below from the extracted repository directory.

## Software

Use R 4.4.1, the reference environment for this package. Core packages are `survival`, `matrixStats`, `mvtnorm`, `statmod`, `MASS`, `jsonlite`, and `digest`. The worker backend uses R's `parallel` package. Recorded versions are in [the reference environment](docs/ACCEPTED_SOFTWARE.json) and [the package version list](docs/VALIDATION_ENVIRONMENT.csv).

Install the core packages in R if needed:

```r
install.packages(c("survival", "matrixStats", "mvtnorm", "statmod", "MASS", "jsonlite", "digest"))
```

Optional R-MAP analyses need JAGS and the R packages `rjags` and `coda`. The additional ALS execution profile also uses `foreach`. No packages are installed automatically. See [optional inputs and dependencies](external/README.md) before running these analyses.

## Small example

```sh
Rscript --vanilla scripts/run_example.R
Rscript --vanilla scripts/summarize_results.R --input results/example --out results/example-summary
```

The example uses one main simulation setting with reduced computation to demonstrate execution. It writes analysis output to `results/example/` and a summary CSV and PDF figure to `results/example-summary/`. Its computation settings are in [examples/smoke.json](examples/smoke.json); [example notes](examples/expected_behavior.md) explain the expected files. It is not a full-study replication.

Use `--out results/my-example` with `run_example.R` to choose another output location. Output directories must be new.

## Choose a scenario

The [scenario guide](docs/SCENARIO_COVERAGE.csv) maps the simulation families to configurations, entry points, and manuscript sections.

| Family | Settings | Configuration | Manuscript |
| --- | --- | --- | --- |
| Main | Same/Diff; N=200 or 300; r=1, 2, or 3; Delta=0, -0.25, or -0.5 | [Main scenarios](config/main/scenarios.json) | Sections 4.1–4.3; Appendices C–D |
| H1–H5 | External baseline-hazard mismatch, non-Weibull baseline, delayed treatment effect, weaker prognostic signal, and heavier censoring | [Robustness configurations](config/H1-H5/) | Sections 4.1 and 4.4; Appendix E |
| A–D | External sample size/coefficient magnitude, covariate correlation, source-coefficient perturbation, and covariate shift | [Sensitivity cases](config/A-D/cases.csv) | Sections 4.1 and 4.4; Appendix E |

List the available scenario IDs:

```sh
Rscript --vanilla scripts/run_simulation.R --list
```

For example, `Same_N300_r2_D0` selects Same external-source structure, trial N=300, a 2:1 treatment-to-control allocation, and Delta=0. The manuscript's 36-setting main study uses the median RECaST-Surv analysis as primary. Detailed median/mean comparisons use 18 selected settings; focused weighting and calibration studies use further subsets. The runner exposes additional outputs, so an available method/scenario combination should not be read as a claim that it appears in every manuscript table. See [analyses and outputs](docs/METHODS_AND_LIMITATIONS.md) for the correspondence.

## Run a scenario or a roster

Run one replication using the full computation settings:

```sh
Rscript --vanilla scripts/run_simulation.R --scenario Same_N300_r2_D0 --reps 1 --seed 902106 --workers 1 --out results/main-one
```

| Option | Purpose |
| --- | --- |
| `--scenario ID` | Select one scenario from `--list`. |
| `--family main`, `H1-H5`, `A-D`, or `all` | Select a family instead of a single scenario. |
| `--reps N` | Set the replication count. |
| `--seed N` | Set a positive integer RNG seed for the new run. |
| `--workers N` | Set the number of workers. |
| `--out PATH` | Write results to a new directory. |
| `--dry-run` | Display the run plan without generating data or fitting models. |
| `--generate-only` | Save generated data without running the analyses. |

[Full-study defaults](config/study_defaults.json) are separate from the reduced example settings. To inspect the complete configured roster at 1,000 replications per setting:

```sh
Rscript --vanilla scripts/run_simulation.R --family all --reps 1000 --seed 902107 --workers 4 --out results/full-roster --dry-run
```

Remove `--dry-run` to execute the plan after choosing the replication count and workers for your resources. The main/H overlap described in Appendix E is represented by shared scenario IDs in the plan. Full rosters require substantial computation; Section 4.2 reports approximately 57.4 seconds per trial for the paper's original B=100 routine under its computing conditions. Runtime on another system can differ.

## Summarize new results

```sh
Rscript --vanilla scripts/summarize_results.R --input results/main-one --out results/main-one-summary
```

Each run directory contains a run manifest and per-replication subdirectories. Main simulation output includes `inputs.rds`, `analysis.rds`, `cox.rds`, `trial_results.csv`, and `status.rds`. The summary command writes `summary.csv`, a combined `trial_results.csv`, `rejection_rates.pdf`, and `SUMMARY_NOTES.txt`. These describe your new run; they are not copies of manuscript tables.

The package also provides [supplementary analysis entries](docs/METHODS_AND_LIMITATIONS.md#supplementary-analysis-entries) for weighting, calibration, and source-coefficient summaries.

## Optional analyses

R-MAP/C-RMAP and ALS have additional requirements listed in [optional inputs and dependencies](external/README.md). The [ALS guide](modules/ALS/README.md) documents the synthetic-data diagnostic entry associated with Section 5 and Appendix F. The archive alone is not sufficient to recreate every paper table.

New simulation estimates may differ numerically from the paper because of random sampling, finite computation, and software environments. For statistical assumptions, interpretation, and numerical considerations, see Sections 3 and 6 and Appendices A–C of the manuscript. [Execution notes](docs/VALIDATION.md) describe the scope of the available checks.

## Citation

Please cite “RECaST-Surv: A Calibrated Borrowing Method for Survival Endpoints in Unequal Randomized Trials,” manuscript version November 24, 2026. Author metadata is provided in [CITATION.cff](CITATION.cff). No final publication DOI is asserted in this package.

## License

This software is released under the [MIT License](LICENSE). Copyright (c) 2026 Dehua Bi. See [license and attribution information](docs/LICENSE_STATUS.md).
