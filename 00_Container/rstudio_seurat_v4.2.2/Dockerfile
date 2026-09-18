# Based on rocker - https://github.com/rocker-org/rocker-versioned
FROM rocker/tidyverse:4.2.2

MAINTAINER Loréna Quatreville (lorena.quatreville@univ-amu.fr)



# ##############################
# R packages system dependencies
# ##############################

RUN apt-get update && apt-get install -y optipng \
                                         pandoc \
                                         libpng-dev \
                                         zlib1g-dev \
                                         libssl-dev \
                                         libcurl4-openssl-dev \
                                         libxml2-dev


# ##########
# R PACKAGES 
# ##########

#### BiocManager
RUN Rscript -e 'install.packages( "BiocManager")'  # currently version 3.16


#### RMarkdown & Reporting
RUN Rscript -e 'install.packages( "knitr")'
RUN Rscript -e 'install.packages( "rmarkdown")'
RUN Rscript -e 'install.packages( "bookdown")'
RUN Rscript -e 'install.packages( "htmltools")'      # browsable


#### Figures & layout management
# ggplot2
RUN Rscript -e 'install.packages( "ggplot2")'
RUN Rscript -e 'install.packages( "cowplot")'        # plot_grid, themes, ...
RUN Rscript -e 'install.packages( "ggpubr")'         # add_summary, geom_signif, ...
RUN Rscript -e 'install.packages( "ggrepel")'        # geom_text_repel, geom_label_repel
RUN Rscript -e 'install.packages( "gridExtra")'      # grid.arrange, ...
RUN Rscript -e 'install.packages( "hexbin")'         # stat_binhex
RUN Rscript -e 'BiocManager::install( "patchwork")'  # +/ operators for ggplots

# plotly
RUN Rscript -e 'install.packages( "plotly")'

# general
RUN Rscript -e 'install.packages( "gplots")'         # heatmap.2
RUN Rscript -e 'install.packages( "heatmaply")'      # heatmaply (interactive)
RUN Rscript -e 'BiocManager::install( "iheatmapr")'  # iheatmap (interactive, uses plotly), dependencies OK with BiocManager
RUN Rscript -e 'install.packages( "pheatmap")'       # pheatmap
RUN apt-get install -y libxt-dev libproj-dev
RUN Rscript -e 'BiocManager::install( "EnhancedVolcano")'
RUN Rscript -e 'BiocManager::install( "ComplexHeatmap")'

# ComplexeHeatmap
##' - The 'magick' package is recommended for better rasterization (install with `install.packages("magick")`). #v2.7.4
RUN apt-get install -y \
    imagemagick \
    libmagick++-dev
#
RUN Rscript -e 'BiocManager::install(
    "UCell",
    ask = FALSE,
    update = FALSE
)' # for module scoring    
    
#### Reporting
RUN Rscript -e 'install.packages( "DT")'             # datatable
RUN Rscript -e 'install.packages( "htmlwidgets")'    # JS (for datatable)
RUN Rscript -e 'install.packages( "kableExtra")'     # kable_styling
RUN Rscript -e 'install.packages( "pander")'         # pander


#### General
RUN Rscript -e 'install.packages( "digest")'         # digest (hashing)
RUN Rscript -e 'install.packages( "dplyr")'
RUN Rscript -e 'install.packages( "forcats")'        # fct_inorder (factors)
RUN Rscript -e 'install.packages( "funr")'           # get_script_path
RUN Rscript -e 'install.packages( "future")'         # plan (multicore)
RUN Rscript -e 'install.packages( "fs")'             # path_sanitize
RUN Rscript -e 'install.packages( "magrittr")'
RUN Rscript -e 'install.packages( "reshape")'        # melt
RUN Rscript -e 'install.packages( "scales")'         # hue_pal
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'install.packages("remotes")'


#### Technology specific
RUN Rscript -e 'install.packages( "Matrix")'         # Force install Matrix 1.5.3 instead of 1.5.1 that cause a bug in FindNeighbors function
RUN Rscript -e 'remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))'
RUN Rscript -e 'remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))' # the one installed is 4.3.0
RUN Rscript -e 'remotes::install_github("satijalab/seurat-data", ref = "v0.2.1")'
RUN Rscript -e 'remotes::install_github("satijalab/azimuth", ref = "v0.4.4")'
RUN Rscript -e 'remotes::install_github("dputhier/scigenex")'

#### Functional analysis
RUN Rscript -e 'BiocManager::install( "org.Hs.eg.db" )' # Human annotation reference
RUN apt-get update && apt-get install -y libxml2-dev libigraph*
RUN Rscript -e 'BiocManager::install( "clusterProfiler")'
RUN Rscript -e 'BiocManager::install( "rrvgo")'
RUN Rscript -e 'BiocManager::install("biomaRt")'
RUN Rscript -e 'BiocManager::install("GenomeInfoDbData")'
RUN Rscript -e 'BiocManager::install("SingleR")'
RUN Rscript -e 'BiocManager::install("celldex")'
RUN Rscript -e 'BiocManager::install("scRNAseq")'
RUN Rscript -e 'BiocManager::install("enrichplot")'
RUN Rscript -e 'BiocManager::install("ReactomePA")'
RUN Rscript -e 'install.packages("gProfileR")'
RUN Rscript -e 'install.packages("gprofiler2")'

# Network analysis
# STRINGdb for protein-protein interaction networks
RUN Rscript -e 'BiocManager::install("STRINGdb", version = "2.10.1")' # 2.10.1

#### Custom dependencies

RUN apt-get update && apt-get install -y curl \
                                         libbz2-dev \
                                         libjpeg-dev \
                                         liblzma-dev \
                                         libxml2-dev

RUN Rscript -e 'BiocManager::install( "GenomicRanges")'
RUN Rscript -e 'BiocManager::install( "rtracklayer")'
RUN Rscript -e 'BiocManager::install( "topGO")'
RUN Rscript -e 'BiocManager::install( "Gviz")'
RUN Rscript -e 'BiocManager::install("MAST")' # FindMarkers statistical test (instead of Wilcox) to account for differences in cell proportions between samples to be compared
RUN Rscript -e 'BiocManager::install("rlist")' # to handle multiple lists

# Additional manual investigation tools (pcaExplorer) 
RUN Rscript -e 'BiocManager::install( "pcaExplorer")' # Also depends on libigraph installed earlier for clusterProfiler  

RUN apt-get update && apt-get install -y libhdf5-dev \
                                         patch
RUN Rscript -e 'install.packages("hdf5r")'
RUN Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'

# Add doubletFinder
RUN Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'

# Add dependencie SciGeneX : mcl
RUN apt-get install -y bzip2
ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"
RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
RUN /bin/bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    echo "export PATH=/opt/conda/bin:$PATH" > /etc/profile.d/conda.sh
ENV PATH /opt/conda/bin:$PATH

RUN conda install -y bioconda::mcl

RUN Rscript -e 'remotes::install_github("satijalab/azimuth", ref = "v0.4.6")'
RUN Rscript -e 'install.packages("scCustomize")'
RUN Rscript -e 'remotes::install_github("dputhier/scigenex@v1.6.1")' # version 1.4.11 is the one installed (for Seurat v4 availability)

RUN apt-get update && apt-get install -y python3 python3-pip  #v3.12.4 #pip 24.0 from /opt/conda/lib/python3.12/site-packages/pip (python 3.12)

RUN pip3 install numpy 
RUN pip3 install pandas
RUN pip3 install leidenalg # v0.4.3
RUN Rscript -e 'install.packages("dbplyr")'
RUN R -e 'devtools::install_github("RobinHankin/Brobdingnag")'

####### Added LATER ON ##########


 # Dependencies for DoubletDecon
RUN Rscript -e 'BiocManager::install("DeconRNASeq", version = "3.16")'
# Dependencies for ShinyFiles
RUN Rscript -e 'remotes::install_version("shiny", version = "1.7.4", repos = "https://cran.r-project.org")'
RUN Rscript -e 'remotes::install_version("fs", version = "1.6.2", repos = "https://cran.r-project.org")'
# RUN Rscript -e 'remotes::install_version("tibble", version = "3.2.1", repos = "https://cran.r-project.org")' # optional if not already compatible
RUN Rscript -e 'remotes::install_version("shinyFiles", version = "0.9.3", repos = "https://cran.r-project.org")'


RUN Rscript -e 'remotes::install_version("R.utils", version = "2.13.0", repos = "https://cran.r-project.org")'
RUN Rscript -e 'remotes::install_version("shinycssloaders", version = "1.1.0", repos = "https://cran.r-project.org")'
RUN Rscript -e 'remotes::install_version("V8", version = "6.0.6", repos = "https://cran.r-project.org")'

# Dependencies for MCL
RUN Rscript -e 'remotes::install_version("expm", version = "1.0-0", repos = "https://cran.r-project.org")'
RUN Rscript -e 'remotes::install_version("MCL", version = "1.0", repos = "https://cran.r-project.org")'

RUN Rscript -e 'BiocManager::install("mygene", version = "3.16")'

RUN Rscript -e 'remotes::install_version("shiny", version = "1.11.1", repos = "https://cran.r-project.org")'
RUN Rscript -e 'devtools::install_github("wleepang/shiny-directory-input@7033725b7737e8c8bd94751a016efb227252c32e")'

RUN Rscript -e 'remotes::install_github("EDePasquale/DoubletDecon@aeee66b", dependencies = FALSE)'



# ICELLNET v2  
# Install ICELLNET dependencies
RUN Rscript -e 'install.packages(c("jetset", "psych"))' \ # v3.4.0 and v2.2.9
RUN Rscript -e 'BiocManager::install("hgu133plus2.db", ask = FALSE, update = FALSE)'
RUN Rscript -e 'remotes::install_github("soumelis-lab/ICELLNET@e10ee4ae5feaf8130764dcb0736957d2a90b828c", subdir="icellnet", dependencies = FALSE)'

# RMCPlus (for Mann-Whitney U test)
# PMCMRplus dependencies (pinned Debian 11 versions)
RUN apt-get update && apt-get install -y \
    libgmp-dev=2:6.2.1+dfsg-3ubuntu1
    libmpfr-dev=4.1.0-3build3
RUN Rscript -e 'remotes::install_version("gmp", version = "0.6-2", repos = "https://cran.r-project.org")'
RUN Rscript -e 'remotes::install_version("Rmpfr", version = "0.8-3", repos = "https://cran.r-project.org")'
# Install PMCMRplus pinned version (for example 1.9.12)
RUN Rscript -e 'remotes::install_version("PMCMRplus", version = "1.9.12", repos = "https://cran.r-project.org")'

# DESeq2 + ZINB-WaVE ecosystem
RUN Rscript -e 'BiocManager::install(c(
    "DESeq2",
    "zinbwave",
    "SingleCellExperiment",
    "SummarizedExperiment",
    "BiocParallel"
    ), ask = FALSE, update = FALSE)'
    
    
RUN Rscript -e 'BiocManager::install("MAST", ask = FALSE, update = FALSE)'

# To have Arial Font
RUN mkdir -p /usr/local/share/fonts/arial
COPY fonts/*.ttf /usr/local/share/fonts/arial/
RUN fc-cache -f -v

# For BetaBinomial model
RUN R -e 'install.packages("remotes", repos="https://cloud.r-project.org"); \
          remotes::install_version("glmmTMB", version="1.1.13", repos="https://cloud.r-project.org")'
