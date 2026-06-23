# IVETTE

IVETTE is the public code repository for the proteomics and clinical analysis pipeline used in the IVETTE COVID-19 study (manuscript in preparation).  
It contains the R scripts and environment definitions needed to reproduce the main analyses, figures, and validation workflows described in the paper.

> **Repository:** [rheumaresearch/IVETTE](https://github.com/rheumaresearch/IVETTE)  
> **Focus:** DIA proteomics discovery, validation experiments and concordance analyses for COVID-related research.

---

## Repository structure

The top-level layout of the project is:

```text
IVETTE/
├── Dockerfile          # Docker image for a reproducible R-based analysis environment
├── LICENSE             # License governing use and redistribution of the code
├── ivette.yml          # Conda environment specification for R and dependencies
└── code/               # Main analysis scripts and modular sub-workflows
    ├── 0.analysis_main.R
    ├── A.Read_data/
    ├── G.Discovery_DIA_FP_MaxLFQ/
    ├── M.Validation_exp/
    ├── O.Concordance_DIA_PRM_complete/
    ├── P.Multiv_predictor/
    ├── S.ELISA_validation/
    └── functions/
```

At a high level:

- **Dockerfile** – Defines a containerized R environment (based on `r-base`) with all required CRAN/Bioconductor packages preinstalled, ensuring that the IVETTE pipeline can be run reproducibly across systems.
- **LICENSE** – Specifies the legal terms under which the code can be used, shared, and modified.
- **ivette.yml** – Provides a Conda environment specification for the IVETTE analysis stack, including R and key dependencies; useful for users who prefer Conda instead of Docker.
- **code/** – Contains the main R script orchestrating the analysis and separate subdirectories for each major analysis component, plus shared utility functions.

---

## Analysis workflow overview

### Main driver script

- **`code/0.analysis_main.R`** – Central “driver” script that sources functions, orchestrates the individual modules (data reading, discovery analysis, validation, concordance, multivariable modeling, ELISA validation), and controls the overall flow of the COVID study pipeline.

### Data import

- **`code/A.Read_data/`** – Scripts related to reading and preparing raw and preprocessed data (e.g. proteomics outputs, clinical variables, metadata tables).  
  This module is responsible for harmonizing input formats and producing clean, analysis-ready data structures for downstream steps.

### Discovery DIA proteomics (FP MaxLFQ)

- **`code/G.Discovery_DIA_FP_MaxLFQ/`** – Discovery-phase DIA proteomics analysis using FP / MaxLFQ-style quantification.
  Typical tasks include quality control, normalization, differential abundance testing, and generation of candidate protein panels for subsequent validation.

### Validation experiments

- **`code/M.Validation_exp/`** – Scripts for validation experiments, such as targeted proteomics runs, replication cohorts, or experimental follow-up designed to confirm findings from the discovery analysis.

### Concordance between DIA and PRM

- **`code/O.Concordance_DIA_PRM_complete/`** – Concordance analysis between discovery DIA data and targeted PRM (parallel reaction monitoring) measurements.


### Multivariable predictor

- **`code/P.Multiv_predictor/`** – Multivariable modeling component. This module generates discriminative models and performance metrics used in the COVID paper.

### ELISA validation

- **`code/S.ELISA_validation/`** – ELISA-based validation scripts, handling data import, normalization, basic statistics, and visualization for enzyme-linked immunosorbent assay experiments.

---

### Shared function library

The **`code/functions/`** directory contains utility scripts that implement reusable functions used across modules. These utilities are typically sourced at the beginning of `0.analysis_main.R` and within module scripts to avoid duplication and keep the pipeline modular.

---

## Environment and reproducibility

IVETTE is designed to be run in a controlled R environment to guarantee reproducibility:

- The **Dockerfile** builds an image encapsulating R, Bioconductor, and CRAN packages used in the project, so users can run the pipeline in a container without configuring the system manually.
- The **`ivette.yml`** file defines a Conda environment with `r-base` and the relevant R packages, providing an alternative for users who prefer Conda-based setups.

Together, these environment definitions ensure that the code in IVETTE faithfully reproduces results regardless of platform, as long as input data are prepared according to the study’s specifications.

---

## How to use this repository

1. **Clone the repository**

   ```bash
   git clone https://github.com/rheumaresearch/IVETTE.git
   cd IVETTE
   ```

2. **Set up the environment**

   - Build the Docker image from `Dockerfile`, or  
   - Create the Conda environment from `ivette.yml`.

3. **Prepare input data**

   - Organize discovery DIA, PRM/validation, ELISA, and clinical data according to the study’s specifications.
   - Update any paths or configuration options referenced by `0.analysis_main.R` and the module scripts.

4. **Run the main script**

   ```bash
   Rscript code/0.analysis_main.R
   ```

   This orchestrates the modules under `code/` to replicate the analyses described in the Rheuma COVID-19 paper, producing results and figures in the configured output locations.

---

## Contact

For questions about IVETTE or its role in the COVID study, please contact the Rheuma research team via the GitHub repository’s issue tracker or project maintainers as listed in the repository metadata.