# ALS registry-informed trial emulation

Section 5 of the manuscript describes an ALS registry-informed emulation using 2,790 synthetic patient records. It considers random and biased trial-cohort sampling and three outcome scenarios: original null, hybrid Delta=0 diagnostic, and alternative Delta=-0.5. Appendix F reports covariate balance and estimated effect ratios.

This guide documents the `modern` profile for the hybrid Delta=0 diagnostic and its baseline/covariate-balance output. It does not identify this profile as a single workflow for regenerating all reported ALS results. Input requirements are specified in [optional inputs and dependencies](../../external/README.md).

## Run the diagnostic

The entry uses `DIAGNOSTIC` for the hybrid Delta=0 scenario, and `random` or `biased` for trial-cohort sampling. To inspect a one-replication plan:

```sh
Rscript --vanilla scripts/run_als.R --profile modern --sampling random --scenario DIAGNOSTIC --input external/data_rwe.csv --reps 1 --seed 902110 --out results/als-diagnostic --dry-run
```

Run this command from the repository root. After supplying the required input, remove `--dry-run` to execute it. Replace `--sampling random` with `--sampling biased` for the other sampling setting, and choose a new output directory. `--reps` and `--seed` control the new run; this entry runs sequentially.

The configurations are [modern/config.R](modern/config.R) and [modern/smd_config.R](modern/smd_config.R). Results are saved as `rep_0001/analysis.rds` and corresponding replication files. Optional `--model-file PATH` adds R-MAP output and a `scores.csv` file when its dependency is supplied. Interpret the hybrid diagnostic and balance output using Section 5.1 and Appendix F.1 of the manuscript.
