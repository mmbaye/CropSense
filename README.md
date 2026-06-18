# CropSense <img src="man/figures/logo.png" align="right" height="80"/>

<!-- badges -->
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**CropSense** is an R package for Near-Infrared (NIR) spectroscopy analysis and
nutritional trait prediction in crop breeding programs. It combines classical
chemometrics (PLS-R) with Gaussian Process uncertainty quantification to tell
you not only *what* a sample contains, but *how confident* the prediction is —
and whether you should send it back to the wet lab.

---

## Supported crops & traits

| Crop | Traits |
|------|--------|
| Sorghum | Protein DW (%), Starch DW (%) |
| Groundnut grain | Protein DW (%), Oil DW (%), Moisture (%) |
| Groundnut flour | Protein DW (%), Oil DW (%), Moisture (%) |

---

## Installation

```r
# From GitHub (recommended)
# install.packages("remotes")
remotes::install_github("your_name/CropSense")
```

---

## Quick start

### Option 1 — Interactive Shiny app

```r
library(CropSense)
run_cropsense()
```

### Option 2 — Command line

```r
library(CropSense)

# Predict traits from a CSV file
results <- predict_traits("my_spectra.csv")
head(results)
```

---

## Input file format

Your CSV or Excel file must contain **exactly these columns**:

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | character | Unique sample identifier |
| `genotype` | character | Genotype / variety name |
| `croptype` | character | Crop type (see supported values below) |
| `wl_400` | numeric | Absorbance at 400 nm |
| `wl_400.5` | numeric | Absorbance at 400.5 nm |
| … | … | … |
| `wl_2499.5` | numeric | Absorbance at 2499.5 nm |

**Supported `croptype` values:**

```
sorghum, sorgho                          → Sorghum model
groundnut, groundnut_grain, peanut, grain → Groundnut grain model
groundnut_flour, flour, farine            → Groundnut flour model
```

**Tip:** Validate your file before predicting:

```r
result <- validate_input("my_spectra.csv", stop_on_error = FALSE)
print_validation(result)
```

---

## Output

```
sample_id | genotype | croptype | protein_pred | proteine_CI95_lower | ... | wetlab_alert
S001      | Tx2741   | sorghum  | 11.23        | 10.45               | ... | FALSE
S002      | BM30     | sorghum  | 12.87        | 9.12                | ... | TRUE   ← send to wet lab
```

Rows where `wetlab_alert = TRUE` have a 95% confidence interval wider than
your threshold — these samples should be confirmed by wet lab analysis.

---

## Architecture

```
CropSense/
├── R/
│   ├── validate_input.R    ← input checking + error messages
│   ├── preprocess.R        ← SNV, MSC, Detrend, SG derivatives
│   ├── predict_traits.R    ← main predict_traits() function
│   ├── gp_regression.R     ← GP LOO + optimization + prediction
│   ├── build_models.R      ← calibration model builder (maintainers)
│   └── run_app.R           ← Shiny app launcher
├── inst/
│   ├── app/app.R           ← full Shiny application (CropSense v2)
│   └── extdata/models/     ← bundled PLS-R calibration models (.rds)
└── tests/testthat/         ← unit tests
```

---

## Citation

If you use CropSense in a publication, please cite:

```
Mbaye, M. (2025). CropSense: NIR Spectroscopy Trait Prediction for Crops.
R package version 0.1.0. https://github.com/your_name/CropSense
```

---

## License

MIT © Modou Mbaye
