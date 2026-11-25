# Analyses and outputs

The associated manuscript, version November 24, 2026, defines RECaST-Surv in Sections 3.2–3.4. The pointwise-median survival summary is primary; the predictive mean is a supplementary analysis. Cox-U denotes the treatment-only concurrent-trial Cox model, and Cox-A denotes the covariate-adjusted model. R-MAP and C-RMAP are the borrowing comparators described in Section 4.1. See the manuscript for their statistical definitions and interpretation.

## Simulation families and paper coverage

The main design comprises 36 combinations of external-source condition, sample size, allocation ratio, and treatment effect (Section 4.1). The main runner includes median and mean RECaST-Surv output, with calibrated and uncalibrated testing labels, and Cox-U/Cox-A output. The manuscript's detailed comparison of the survival summaries concerns 18 selected settings (Section 4.2 and Appendix C). The availability of additional code output does not expand the set of results reported in the paper.

Robustness scenarios H1–H5 examine external baseline hazards, a piecewise-exponential baseline, delayed treatment effect, weaker prognostic coefficients, and heavier potential dropout. The paper reports 14 configurations with the shared trials described in Table A12. Tables A12–A14 present testing, point-estimation, and median confidence-interval results. Mean-summary output is also available in the runner, but is not presented as a reported H1–H5 comparison here.

Sensitivity Sets A–D contain 29 designs, each with Delta=0 and -0.5. Their configurations retain the individual external-size, correlation, coefficient-perturbation, and support-restriction cases. The RECaST-Surv branch for these sets uses the median summary. The paper reports testing results in Tables A15–A16. The [scenario guide](SCENARIO_COVERAGE.csv) identifies the relevant files and sections.

## Reading the output

| Output label | Meaning |
| --- | --- |
| `scenario_id`, `replication_id`, `method` | Identify the setting, replication, and analysis. |
| `RECaST_median_*`, `RECaST_mean_*` | Identify the pointwise survival summary and calibrated or uncalibrated testing output. |
| `estimate` | The method's estimated effect ratio; RECaST-Surv uses the O/E summary in Section 3.3. |
| `p`, `cutoff`, `reject` | The reported test quantities and decision. |
| `CI_status`, `CI_lower`, `CI_upper`, `CI_length` | Availability and endpoints/length of the reported interval. |
| `reference`, `reference_id` | The evaluation reference where one is supplied for that method and setting. |

`summary.csv` includes available counts, rejection proportions and Monte Carlo standard errors, point summaries, and interval summaries where applicable. Calibrated and uncalibrated RECaST-Surv testing labels share the corresponding point estimate and interval, as described in Section 4.1. The paper uses method-specific evaluation references; the output's `NA` values should be retained when a reference or interval is unavailable. The delayed-effect H3 setting and its reference-based summaries are discussed in Section 4.1 and Appendix E.

Use the manuscript's distinctions among the O/E summary, Cox estimates, and R-MAP estimates when interpreting output. The software's output columns do not imply that all methods estimate a common parameter or that every column is reported for every scenario.

## Supplementary analysis entries

These entries operate on saved output from new full-computation runs. Their focused simulation settings are documented in [config/focused/approved.json](../config/focused/approved.json).

| Entry | Use | Manuscript |
| --- | --- | --- |
| `scripts/run_sensitivity.R --kind weights` | Original, equal, and normalized source weighting. | Section 3.4; Appendix C.4, Table A6 |
| `scripts/run_sensitivity.R --kind cv5` | Five-fold cross-validated source weighting. | Section 3.4; Appendix C.4, Table A6 |
| `scripts/run_sensitivity.R --kind calibration` | Repeated-splitting B=100/200/500 and the full-size bootstrap comparison. | Section 4.2; Appendix C.3, Table A5 |
| `scripts/report_focused.R --kind full-vs-split` | Full-control and split-control score-distribution summaries. | Appendix C.3, Figure A1 and Table A4 |
| `scripts/report_focused.R --kind coefficients` | Source-coefficient error and association summaries. | Section 4.2; Appendix C.4 |

For example, plan a weighting analysis after generating the corresponding full-computation trial:

```sh
Rscript --vanilla scripts/run_sensitivity.R --analysis results/main-one/Same_N300_r2_D0/rep_0001/analysis.rds --kind weights --seed 902112 --out results/weighting --dry-run
```

`run_sensitivity.R` requires `--analysis`, `--kind`, `--seed`, and `--out`. The weighting entries use Same, N=300, r=2 at Delta=0 or -0.5; the calibration entry uses the null setting. They do not accept reduced-example output. `report_focused.R` requires `--kind`, `--input`, and `--out`; its complete-bank summaries require 1,000 null replications and, for coefficient summaries, 1,000 alternatives. Use a new output directory for each analysis. The entries support `--dry-run` for planning.

For statistical assumptions, scope, and numerical limitations, consult Sections 3 and 6, Appendix A, and Appendix C.5. This guide does not replace those discussions.
