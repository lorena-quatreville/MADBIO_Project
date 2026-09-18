
# R 4.6.x / Bioconductor 3.23
FROM rocker/tidyverse:4.6.0

# Utilisateur RStudio
ENV DEFAULT_USER=quatreville
ENV HOME=/home/quatreville

MAINTAINER Loréna Quatreville (lorena.quatreville@univ-amu.fr)

# ============================================================
# System dependencies
# ============================================================

RUN apt-get update && apt-get install -y --no-install-recommends \
        optipng \
        pandoc \
        libpng-dev \
        zlib1g-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libxml2-dev \
        wget \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*
    
    
# ============================================================
# Bioconductor
# ============================================================

RUN Rscript -e 'install.packages("BiocManager", repos="https://cloud.r-project.org")' \
    && Rscript -e 'BiocManager::install(version = "3.23", ask = FALSE)'


# ============================================================
# RMarkdown / reporting
# ============================================================

RUN Rscript -e 'install.packages(c( \
    "knitr", \
    "rmarkdown", \
    "bookdown", \
    "htmltools" \
), repos="https://cloud.r-project.org")'


# ============================================================
# Figures / plotting
# ============================================================

RUN Rscript -e 'install.packages(c( \
    "ggplot2", \
    "cowplot", \
    "ggpubr", \
    "ggrepel", \
    "gridExtra", \
    "hexbin", \
    "plotly", \
    "DT", \
    "htmlwidgets", \
    "kableExtra", \
    "pander" \
), repos="https://cloud.r-project.org")'

RUN Rscript -e 'BiocManager::install("patchwork", ask=FALSE)'


# ============================================================
# General R packages
# ============================================================

RUN Rscript -e 'install.packages(c( \
    "digest", \
    "dplyr", \
    "forcats", \
    "funr", \
    "future", \
    "fs", \
    "magrittr", \
    "reshape", \
    "scales", \
    "devtools", \
    "remotes" \
), repos="https://cloud.r-project.org")'


# ============================================================
# Conda / mcl
# ============================================================

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
        -O /tmp/miniconda.sh \
    && /bin/bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

ENV PATH="/opt/conda/bin:${PATH}"

RUN conda tos accept --override-channels \
        --channel https://repo.anaconda.com/pkgs/main \
    && conda tos accept --override-channels \
        --channel https://repo.anaconda.com/pkgs/r

RUN conda install -y --channel conda-forge --channel bioconda --strict-channel-priority mcl \
    && conda clean -afy


# ============================================================
# CmdStanR
# ============================================================

# Current CmdStanR release
RUN Rscript -e 'install.packages( \
    "cmdstanr", \
    repos=c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"), \
    dependencies=TRUE \
)'


# ============================================================
# CmdStan
# ============================================================

RUN Rscript -e 'cmdstanr::check_cmdstan_toolchain(fix=TRUE)' \
    && Rscript -e 'cmdstanr::install_cmdstan(cores=2)'


# ============================================================
# scComp
# ============================================================

RUN Rscript -e 'BiocManager::install("sccomp", ask=FALSE, update=TRUE)'


# ============================================================
# Sanity checks
# ============================================================

RUN Rscript -e ' \
    cat("R version: ", R.version.string, "\n"); \
    cat("Bioconductor: ", as.character(BiocManager::version()), "\n"); \
    cat("sccomp: ", as.character(packageVersion("sccomp")), "\n"); \
    cat("cmdstanr: ", as.character(packageVersion("cmdstanr")), "\n"); \
    cat("CmdStan: ", as.character(cmdstanr::cmdstan_version()), "\n"); \
'
