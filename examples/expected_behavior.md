# Reduced example

From the repository root:

```sh
Rscript --vanilla scripts/run_example.R
```

The example selects Same, N=200, r=1, Delta=0 and uses the reduced settings in [smoke.json](smoke.json). It creates one replication with median/mean RECaST-Surv output and Cox-U/Cox-A output in `results/example/`. The RECaST-Surv interval fields are unavailable because the example uses only two calibration splits; it is intended to demonstrate execution rather than reproduce a paper estimate.

Use [summarize_results.R](../scripts/summarize_results.R) as shown in the main [README](../README.md) to create a CSV summary and PDF figure. Generated results are created locally and are not included in the archive. Full-study settings are in [study_defaults.json](../config/study_defaults.json).
