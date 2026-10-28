# CropSense

**An integrated R Shiny platform for multi-instrument NIR-based quality phenotyping in crop breeding**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/mmbaye/CropSense/blob/main/LICENSE)
[![R](https://img.shields.io/badge/R-≥4.2-blue.svg)]()
[![Version](https://img.shields.io/badge/version-v6-green.svg)]()
[![SoftwareX](https://img.shields.io/badge/Published-in-SoftwareX-brightgreen.svg)]()

---

## What is CropSense?

CropSense is an **open-source, browser-based application** written in R (Shiny) that unifies the full near-infrared (NIR) spectroscopy workflow for crop quality phenotyping:

> **Acquire → Harmonize → Explore → Calibrate & Predict → Report**

It was developed at the **Senegalese Agricultural Research Institute (ISRA/CERAAS)**, Dakar, to support routine nutritional quality phenotyping in breeding programs for **groundnut, sorghum, pearl millet, cowpea, fonio, and sesame**.

CropSense is not limited to agriculture: any laboratory operating multi-instrument NIR spectroscopy (food, pharmaceutical, environmental, botanical) can adapt the platform.

### Key capabilities

| Stage | What CropSense does |
|-------|-------------------|
| **Acquire** | Imports spectral files from **FOSS**, **Bruker FT-NIR**, and **MicroNIR** spectrometers (CSV, TSV, TXT, XLS, XLSX) |
| **Harmonize** | Auto-detects delimiters, normalises wavelengths, extracts metadata, maps sample–reference links |
| **Explore** | Interactive spectral overlays (Plotly), preprocessing (SNV, MSC, detrending, Savitzky–Golay, 1st/2nd derivatives), PCA score & loading plots |
| **Calibrate** | PLS-R calibration with cross-validation; Gaussian process regression with 95 % predictive uncertainty |
| **Predict** | Compatibility-checked application of registered models to new spectra (instrument, signal, wavelength grid, completeness) |
| **Report** | 5 automated report types in self-contained HTML (printable to PDF) + CSV exports |

### The Model Registry

A central feature of CropSense v6 is the **persistent Model Registry**:

- Fitted calibrations (PLS-R or GP) are stored with their **exact preprocessing configuration**, instrument metadata, wavelength grid, and validation diagnostics.
- **Grouped biological-sample validation**: technical replicates of the same biological sample are never split across train/holdout partitions.
- **Status lifecycle**: Candidate → Active → Archived.
- Backed by **SQLite** (metadata) + **RDS** (model objects). Persists across sessions and sampling campaigns.

---

## Installation

### Prerequisites

- **R** ≥ 4.2: [https://www.r-project.org/](https://www.r-project.org/)
- A modern web browser (Chrome, Firefox, Edge, or Safari)
- No internet connection required at runtime (local Shiny server)

### Install dependencies

```r
# Core packages
if (!require("remotes")) install.packages("remotes")

remotes::install_cran(c(
  "shiny", "bslib", "plotly", "ggplot2",
  "pls", "prospectr", "MASS",
  "readr", "readxl", "dplyr", "tidyr",
  "DBI", "RSQLite", "digest", "testthat", "lintr"
))
```

### Clone and run

```bash
git clone https://github.com/mmbaye/CropSense.git
cd CropSense
```

```r
# In R or RStudio
setwd("/path/to/CropSense")
shiny::runApp("app.R")
```

CropSense will open automatically in your default browser at `http://localhost:3838` (or the port shown in the R console).

> **Tip:** If you work on Windows and encounter path issues, use forward slashes (`/`) instead of backslashes (`\`) in `setwd()`.

---

## Quick start

### 1. Import a spectral file

1. Open the **Import** module.
2. Select a CSV, TSV, XLS, or XLSX file containing your spectral data (samples in rows, wavelengths in columns).
3. CropSense automatically:
   - Detects the file delimiter
   - Identifies and normalises wavelength columns
   - Preserves non-spectral metadata (genotype, batch, reference traits)
4. Review the dataset summary (n spectra, n channels, wavelength range).

### 2. Explore and preprocess

1. Open the **Spectra** module.
2. Overlay a subset of spectra to check for artefacts.
3. Select a preprocessing method (e.g., SNV → 1st SG derivative, window 11, polynomial 2).
4. Compare preprocessed spectra with the raw curves.

### 3. Run PCA *(optional)*

1. Open the **PCA** module.
2. Select the number of components to display.
3. Colour samples by genotype or batch.
4. Inspect the scree plot and score/loadings.

### 4. Calibrate *(requires reference data)*

1. Open the **PLS Calibration** module.
2. Select the target trait (e.g., protein %).
3. Choose preprocessing and maximum latent variables.
4. Review the RMSECV curve and measured-vs-predicted plot.
5. *(Optional)* Run the **GP & Uncertainty** module for a nonlinear model with predictive intervals.

### 5. Register the model

1. Open the **Model Registry** module.
2. Fill in: project, crop, instrument, biological sample grouping key.
3. Click **Register**. The model is saved as a **Candidate**.

### 6. Predict new samples

1. Open the **Predictions** module.
2. Select the registered model from the Registry.
3. Import the new spectral file.
4. CropSense checks instrument, signal, wavelength, and completeness.
5. Predictions are generated with domain-check flags (and uncertainty intervals for GP models).

### 7. Generate reports

1. Open the **Reports** module.
2. Select one or more report types.
3. Download the HTML file (open in a browser, print to PDF) and/or the CSV data tables.

---

## Project structure

```
CropSense/
├── app.R                  # Main Shiny app (UI + server composition)
├── modules/
│   ├── import.R           # Import & harmonization module
│   ├── spectra.R          # Spectral visualization & preprocessing
│   ├── pca.R              # PCA module
│   ├── pls.R              # PLS-R calibration module
│   ├── gp.R               # Gaussian process regression module
│   ├── registry.R         # Model Registry (SQLite + RDS)
│   ├── predictions.R      # Prediction of new samples
│   ├── datatable.R        # Non-spectral metadata inspector
│   └── reports.R          # Automated report generation
├── utils/
│   ├── preprocess.R       # SNV, MSC, SG, derivatives
│   ├── wavelength.R       # Wavelength detection & normalization
│   ├── grouped_split.R    # Biological-sample-aware train/holdout split
│   ├── compatibility.R    # Prediction compatibility checks
│   └── report_gen.R       # HTML/CSV report templates
├── tests/
│   └── testthat/
│       ├── test_preprocess.R
│       ├── test_wavelength.R
│       ├── test_grouped_split.R
│       └── test_compatibility.R
├── registry/              # SQLite DB + RDS model files (created at runtime)
├── LICENSE                # MIT License
├── README.md              # This file
└── renv/                  # Package lockfile (reproducible environment)
```

> **Note:** The `registry/` directory is created automatically on first use. Back it up if you want to preserve registered models across machine changes.

---

## Supported instruments

| Instrument | Type | Typical output | Notes |
|------------|------|---------------|-------|
| FOSS (e.g., XDS, Infratec) | Desktop dispersive NIR | CSV/XLS, samples × wavelengths | Transmittance or reflectance |
| Bruker (e.g., MPT, TENSOR) | Desktop FT-NIR | CSV/XLS | May require DataBridge for row-major exports |
| MicroNIR (e.g., 2.0, 2.5) | Portable / handheld | CSV/TXT | Smaller wavelength range, technical replicates common |

> **Important:** Harmonisation of file formats does **not** imply transferability of calibration coefficients between instrument families. CropSense enforces instrument and wavelength-grid compatibility before prediction.

---

## Preprocessing methods

| Method | Description |
|--------|-------------|
| SNV | Standard Normal Variate |
| MSC | Multiplicative Scatter Correction |
| Detrending | Baseline trend removal |
| SNV + Detrending | Sequential |
| SG smoothing | Savitzky–Golay (adjustable window & polynomial) |
| 1st SG derivative | Slope emphasis |
| 2nd SG derivative | Curvature emphasis (noise-sensitive) |

---

## Reporting

CropSense generates five report types as self-contained HTML documents:

1. **Full Technical Report**
2. **Spectral Inventory & Quality Report**
3. **Spectral Exploration Report**
4. **Calibration Performance Report**
5. **Predictions & Uncertainty Report**

Each report can be:
- Opened directly in a browser
- Printed / saved as PDF
- Accompanied by a CSV data export

---

## Testing

Unit tests are provided for core analytical functions:

```r
# Run all tests
devtools::test()

# Or manually
library(testthat)
test_dir("tests/testthat")
```

Covered: preprocessing transformations, wavelength normalisation, grouped split logic, model serialisation, compatibility checks.

---

## Known limitations

- Independent **external validation** is not yet automated (conducted outside the application).
- **Cross-instrument model transfer** is not supported.
- Complex instrument file formats (Bruker row-major, multi-level ICP-MS) require the companion **DataBridge** workflow (in development).
- GP uncertainty intervals are **model-based**, not externally validated prediction coverage.
- Single-user, local-server architecture (no multi-user concurrent access yet).

---

## Roadmap

- [ ] Native DataBridge integration (Bruker, ICP-MS native formats)
- [ ] Persistent institutional spectral library
- [ ] Independent validation workflow
- [ ] Publication-quality figure export
- [ ] Standardised model governance (Candidate → Validated → Operational → Retired)
- [ ] Interoperability with breeding data management systems
- [ ] Multi-user / networked access
- [ ] CropSenseR package (script-based analysis and figure generation)

---

## Citation

If you use CropSense in your research, please cite:

> Modou MBAYE, khalil Kane, Ndeye Coura Fall. **CropSense: An integrated platform for multi-instrument NIR-based quality phenotyping in crop breeding.** *SoftwareX* (forthcoming). https://doi.org/[DOI]

And the software:

> CropSense. Open-source R Shiny platform for NIR-based quality phenotyping. Version v6. GitHub: https://github.com/mmbaye/CropSense. MIT License.

*(Update the citation block above with the published author list, journal volume/issue, and DOI once the paper is accepted.)*

---

## Contributing

Contributions are welcome. Please:

1. Fork the repository.
2. Create a feature branch (`git checkout -b feature/improvement-name`).
3. Add or update unit tests.
4. Ensure `devtools::test()` and `lintr::lint_package()` pass.
5. Submit a pull request with a clear description.

### Issue reporting

Use the [GitHub Issues](https://github.com/mmbaye/CropSense/issues) tracker for bug reports, feature requests, and questions.

---

## License

This project is distributed under the **MIT License**. See [LICENSE](LICENSE) for details.

---

## Acknowledgements

Developed within the crop improvement programs of the **Senegalese Agricultural Research Institute (ISRA/CERAAS)**, Thies, Senegal.

funding: CIMMYT VACS Processing

---

## Contact

**Modou Mbaye** (corresponding author)
ISRA/CERAAS, Dakar, Senegal
Email: modou.mbaye@isra.sn
