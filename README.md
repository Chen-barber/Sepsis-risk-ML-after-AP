# Sepsis-3 risk calculator

Local Shiny web calculator for the reported GAMBoost model.

Enter the seven consensus features to obtain an individualized Sepsis-3 probability. The app shows KernelSHAP waterfall and force plots and can download a one-page result PDF.

## Manuscript wording

To facilitate practical implementation, we developed a web-based risk calculator (available at: https://github.com/Chen-barber/Sepsis-risk-ML-after-aspiration-pneumonia/tree/web-calculator). This user-friendly application not only generates individualized sepsis risk probabilities based on the seven variables but also provides real-time local interpretability through integrated SHAP waterfall and force plots, supporting transparent clinical decision-making.

## Run

```text
cd shiny
Rscript --vanilla run_app.R
```

Default local URL: `http://127.0.0.1:8080` (override with environment variables `SHINY_HOST` and `SHINY_PORT`).

## Dependencies

Install once:

```r
install.packages(c(
  "shiny", "bslib", "mboost", "kernelshap", "shapviz",
  "ggplot2", "patchwork"
))
```

## Input features

- 24-h blood urea nitrogen (mmol/L)
- 24-h calcium (mmol/L)
- 24-h chloride (mmol/L)
- Mechanical ventilation (0/1)
- Norepinephrine use (0/1)
- 24-h platelet count (x10^9/L)
- 24-h systolic blood pressure (mmHg)

Threshold options: training-set Youden cut-off (default) or fixed 0.50.

## Privacy

`shiny.Rdata` contains the fitted GAMBoost object and an anonymous background matrix of seven clinical features for KernelSHAP. No patient identifiers are included.
