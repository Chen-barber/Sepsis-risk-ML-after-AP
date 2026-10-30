# Sepsis-risk-ML-after-AP

Public code deposit for the reported analysis and the interactive risk calculator.

## Contents

- `split_mimic_clean.R` — stratified 7:3 split of the imputed MIMIC analysis set (seed 341)
- `ml_mimic_clean.R` — seven classifiers with fixed hyperparameters (seed 123)
- `shiny/` — local web calculator (GAMBoost probability + SHAP waterfall/force plots)
- `requirements.txt` / `requirements-R.txt` — software versions
- `LICENSE` — MIT Licence

Patient-level data are not included.

## Web calculator

Branch `web-calculator` (and the `shiny/` folder on this repository):

https://github.com/Chen-barber/Sepsis-risk-ML-after-aspiration-pneumonia/tree/web-calculator

Suggested manuscript text:

> To facilitate practical implementation, we developed a web-based risk calculator (available at: https://github.com/Chen-barber/Sepsis-risk-ML-after-aspiration-pneumonia/tree/web-calculator). This user-friendly application not only generates individualized sepsis risk probabilities based on the seven variables but also provides real-time local interpretability through integrated SHAP waterfall and force plots, supporting transparent clinical decision-making.

Run locally:

```text
cd shiny
Rscript --vanilla run_app.R
```

The console prints a local URL (default `http://127.0.0.1:8080`).

## Data access

MIMIC-IV is obtained from PhysioNet after credentialing. Access MIMIC-IV at https://physionet.org/content/mimiciv/ under the PhysioNet Credentialed Health Data Use Agreement. Do not redistribute MIMIC-derived tables via this repository. The external institutional cohort is available only under local ethics approval and may be requested from the corresponding author.

## Environment

Computing environment: R 4.6.0 and Python 3.14.5. Package versions are listed in `requirements-R.txt` and `requirements.txt`.

Please cite the published article when available and acknowledge PhysioNet and MIMIC-IV.
