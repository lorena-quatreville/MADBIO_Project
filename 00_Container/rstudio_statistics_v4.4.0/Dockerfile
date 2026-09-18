# Use a Rocker image with R 4.4 (Required for Bioconductor 3.20)
FROM rocker/tidyverse:4.4.0

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

#### Install BiocManager and set Bioconductor version to 3.20
RUN Rscript -e 'install.packages("BiocManager", repos="https://cran.r-project.org")' \
    && Rscript -e 'BiocManager::install(version = "3.20")'


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
RUN Rscript -e 'remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))' # version 5.1.0 installed
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
RUN Rscript -e 'BiocManager::install("clusterProfiler")'
RUN Rscript -e 'BiocManager::install("ReactomePA")'
RUN Rscript -e 'BiocManager::install("DCATS")'
RUN Rscript -e 'install.packages("gProfileR")'
RUN Rscript -e 'install.packages("gprofiler2")'

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
RUN Rscript -e 'remotes::install_github("dputhier/scigenex")' # version 1.5.3 is installed

RUN apt-get update && apt-get install -y python3 python3-pip 

RUN pip3 install numpy 
RUN pip3 install pandas
RUN pip3 install leidenalg 
RUN Rscript -e 'install.packages("dbplyr")'
RUN R -e 'devtools::install_github("RobinHankin/Brobdingnag")'


##### Package pour faire tourner scDC (Single cell differential composition analysis)

## Some CRAN packages required by scDC

RUN R -e 'install.packages(c("parallel", "DescTools", "lme4", "reshape2", "ggridges", "lme4", "mice"))'

## Some BioConductor packages required by scDC
RUN R -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install("scran")'

## Some Github packages required by scDC
RUN R -e 'if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools"); devtools::install_github("taiyunkim/scClustBench")'


## Installing scDC 
RUN R -e 'devtools::install_github("SydneyBioX/scDC")'


# Add Monocle for trajectory analysis
RUN apt-get update && apt-get install -y libgdal30=3.4.1+dfsg-1build4 libgdal-dev=3.4.1+dfsg-1build4 # necessary for installing Monocle3
RUN Rscript -e "remotes::install_github('cole-trapnell-lab/monocle3@v1.3.1')" # Installs Monocle version 3

RUN R -e "remotes::install_version('R.utils', version = '2.13.0', repos='https://cran.r-project.org', upgrade='never')"
RUN Rscript -e "remotes::install_github("satijalab/seurat-wrappers")" # enable conversion from Seurat to Monocle's cell_data_set (version 0.4.0)

# Add scComp for compositional analysis for flow cytometry data
RUN R -e 'BiocManager::install("sccomp", ask=FALSE, update=FALSE)' # v1.10.0




remotes::install_version("cmdstanr", version = "0.8.1", repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"), upgrade = "never", dependencies = FLASE)



# Fixed CmdStan version
RUN Rscript -e '\
cmdstanr::install_cmdstan( \
    version = "2.36.0", \
    cores = 2, \
    quiet = FALSE \
)'

# Fixed instantiate version
RUN Rscript -e '\
remotes::install_version( \
    "instantiate", \
    version = "0.2.3", \
    repos = "https://cloud.r-project.org", \
    upgrade = "never", \
    dependencies = FALSE \
)'

# sccomp for Bioconductor 3.20
RUN Rscript -e '\
BiocManager::install( \
    "sccomp", \
    ask = FALSE, \
    update = FALSE \
)'
