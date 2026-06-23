FROM r-base:latest
LABEL maintainer="Diego Fuentes <dfuentesp@tauli.cat>"
LABEL description="Reproducible R environment for the project IVETTE, including Bioconductor and CRAN packages."

# System dependencies commonly needed to compile CRAN/Bioconductor packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libpng-dev \
    libjpeg-dev \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor manager and Bioconductor packages
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org'); \
    BiocManager::install(c('affy', 'limma', 'roastgsa'), ask = FALSE, update = TRUE)"

# Install CRAN packages corresponding to the conda-forge R packages:
RUN R -e "install.packages(c( \
    'cluster',       \
    'hwriter',       \
    'multcomp',      \
    'amap',          \
    'gtools',        \
    'quantreg',      \
    'lme4',          \
    'lmerTest',      \
    'psych',         \
    'ggplot2',       \
    'ggrepel',       \
    'eulerr',        \
    'VennDiagram',   \
    'colorspace',    \
    'ppcor',         \
    'randomForest',  \
    'pROC',          \
    'glmnet',        \
    'minpack.lm',    \
    'drc',           \
    'combinat',      \
    'ggpubr',        \
    'rstatix',       \
    'ggtext',        \
    'openxlsx',      \
    'Hmisc',         \
    'gt',            \
    'emmeans',       \
    'gtsummary',     \
    'patchwork',     \
    'Minirand',      \
    'iq',            \
    'remotes'       \
    ), repos = 'https://cloud.r-project.org')"

# Missing dependencies, sinkr from GitHub (marchtaylor/sinkr)
RUN R -e "remotes::install_github('marchtaylor/sinkr')"

# Default workdir for your R scripts
WORKDIR /usr/src/app

# Copy your R scripts into the image (optional; adjust to your repo layout)
COPY ./code/ ./scripts/

# Example default command (override as needed)
CMD [\"Rscript\", \"scripts/0.analysis_main.R\"]