# Stroke Trial Reanalysis

A statistical reanalysis of the **International Stroke Trial** (IST): 19,435 patients with acute ischaemic stroke, randomised in a 2 by 2 factorial design to aspirin or not, and to heparin (low dose, medium dose) or not. The analysis reproduces the published primary result from the patient level data in R, then adds a covariate adjusted estimate, a benefit and harm comparison by heparin dose, survival models, multiplicity corrected subgroups, missing data sensitivity analyses and a cross validated prognostic model. A SAS program runs the same primary analyses.

**[Dashboard](https://satwik-kotta.github.io/stroke-trial-reanalysis/)** · **[Statistical report](https://satwik-kotta.github.io/stroke-trial-reanalysis/report.html)**

## Findings

| Death or dependency at 6 months | Treated | Control | Difference per 100 (95% CI) | p |
|---|---|---|---|---|
| Aspirin, crude | 62.2% | 63.5% | −1.3 (−2.6 to 0.1) | 0.07 |
| Aspirin, adjusted for baseline prognosis | | | −1.5 (−2.6 to −0.4) | 0.016 |
| Heparin, adjusted | 62.9% | 62.9% | 0.0 (−1.3 to 1.2) | 1.00 |

- **The published result reproduces exactly**: 62.2% against 63.5% (Lancet 1997, 2p = 0.07). Adjusting for baseline prognosis gives 1.5 fewer dead or dependent per 100, close to the published adjusted 1.4.
- **Heparin's benefit is cancelled by bleeding.** Per 1,000 patients in 14 days, medium dose prevents 2 recurrent ischaemic strokes but causes 13 haemorrhagic strokes and 17 major extracranial bleeds. Low dose prevents 8 and causes 3 and 2.
- **No subgroup holds up after correction.** Heparin by stroke type gives p = 0.04 uncorrected and 0.29 after Holm correction across 8 factors.
- **The aspirin result survives** multiple imputation (20 imputations, p = 0.016), per protocol, excluding non strokes and scanned patients only. It loses significance only when every missing outcome counts against aspirin.
- **Prognostic model**: cross validated AUC 0.777, calibration slope 0.99.

## Data checks

- 150 patients have no 6 month outcome: 53 coded 8 or 9 and 97 coded 0, a code the documentation does not define. All treated as missing (99.2% follow up).
- The pilot coded medium dose heparin as H (245 patients), merged with M.
- Atrial fibrillation not recorded for 984 pilot patients, kept as their own level in models.
- 420 patients turned out not to have had a stroke; kept as randomised, excluded in a sensitivity analysis.
- Every baseline standardised mean difference between arms is under 0.025.

## Methods

| Step | File |
|---|---|
| Data checks and analysis dataset | `R/prepare.R` |
| Newcombe intervals, standardised (g computation) risk differences with bootstrap, subgroup interaction tests, Holm, Kaplan Meier, Cox with Schoenfeld test, cross validated AUC and calibration | `R/methods.R` |
| All analyses, including multiple imputation with mice | `analysis/run.R` |
| Primary analyses in SAS (PROC FREQ, LOGISTIC, LIFETEST, PHREG) | `sas/ist_primary.sas` |
| Tests on simulated data and a published worked example | `tests/testthat/` |

## Run it

```bash
bash scripts/get_data.sh            # downloads the IST database and writes data/ist_raw.csv
Rscript analysis/run.R              # writes docs/data/results.json and outputs/*.csv
Rscript -e 'testthat::test_dir("tests/testthat")'
```

R packages: survival, mice, jsonlite, data.table, testthat.

## Data

Sandercock P, Niewada M, Czlonkowska A. The International Stroke Trial database (version 2). University of Edinburgh, 2011. https://doi.org/10.7488/ds/104. Open Data Commons Attribution licence. Described in *Trials* 2011;12:101. Trial results: International Stroke Trial Collaborative Group, *Lancet* 1997;349:1569–81.

## Limits

A 1990s trial, before thrombolysis and thrombectomy were routine. A third of patients had no CT scan before randomisation. The results describe aspirin and heparin as tested then, not current practice.
