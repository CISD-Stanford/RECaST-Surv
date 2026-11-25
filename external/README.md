# Optional inputs and dependencies

| Analysis | Required material | Identity requirement |
| --- | --- | --- |
| R-MAP/C-RMAP | JAGS, the R packages `rjags` and `coda`, and the study-specific `rMAP_Surv_Pois_Fit_ctl.txt` model. | The model must match [RMAP_MODEL_REQUIRED.json](RMAP_MODEL_REQUIRED.json). |
| ALS | The separately supplied synthetic `data_rwe.csv` input. | [run_als.R](../scripts/run_als.R) checks the specified input's SHA-256 before reading it. |

The model text and CSV are excluded from this archive. Obtain authorized copies separately and pass their paths to the relevant entry. JAGS and R packages use their own distributions and licenses; installing them does not supply the study-specific model. The repository's software license status is described in [LICENSE_STATUS.md](../docs/LICENSE_STATUS.md).

## R-MAP and C-RMAP entries

Run these commands from the repository root. For example, the following plans a null-scenario R-MAP run:

```sh
Rscript --vanilla scripts/run_rmap.R --scenario Same_N300_r2_D0 --reps 1000 --seed 902108 --workers 1 --out results/rmap-null --model-file external/rMAP_Surv_Pois_Fit_ctl.txt --dry-run
```

After supplying the model, remove `--dry-run` to execute the plan. The runner also accepts `--family main`, `H1-H5`, `A-D`, or `all`. Replication counts can be at most 1,000 in this entry. Output includes `scores.csv`. To use data already saved by a full-computation primary run, additionally supply `--input-run PATH`; otherwise the entry generates its own inputs.

For an alternative-scenario score bank, use `--scenario Same_N300_r2_Dm05`, a different seed, and `--out results/rmap-alternative`. C-RMAP can then read the two score files:

```sh
Rscript --vanilla scripts/calibrate_rmap.R --null results/rmap-null/scores.csv --null-bank Same_N300_r2_D0 --evaluation results/rmap-alternative/scores.csv --evaluation-bank Same_N300_r2_Dm05 --out results/crmap
```

The selected null bank must contain exactly 1,000 complete scores. Use bank labels from your score files. Calibration output is saved as `calibration.rds` and `decisions.csv`. The comparison is described in Section 4.1 of the manuscript; the [ALS guide](../modules/ALS/README.md) covers the separate synthetic-data diagnostic entry.
