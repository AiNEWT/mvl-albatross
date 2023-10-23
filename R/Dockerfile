ARG R_BASE_VERSION=4.3.1
ARG SHINY_VERSION=1.5.19.995

FROM --platform=linux/amd64 rocker/rstudio:${R_BASE_VERSION}

# LABEL org.opencontainers.image.licenses="GPL-2.0-or-later" \
#       org.opencontainers.image.source="[SOME GITHUB REPOs]" \
#       org.opencontainers.image.vendor="Albatross Analyzer Project" \
#       org.opencontainers.image.authors="Albatross Analyzer <albatross@miningvizlab.com>"

# Update OS and Install Build Essential, Library
ARG PREBUILD_DEPS="software-properties-common gnupg2"
ARG BUILD_DEPS="build-essential apt-utils binutils cmake gfortran git libblas-dev libcairo2-dev libcurl4-openssl-dev libjpeg-dev liblapack-dev libnlopt-dev libpng-dev libturbojpeg libxt-dev libcurl4-openssl-dev libssl-dev libxml2-dev pkg-config sudo"

RUN apt-get update \
    && apt-get install -y --no-install-recommends $PREBUILD_DEPS \
    && apt-get update \
    && apt-get install -y --no-install-recommends $BUILD_DEPS \
    && rm -rf /var/lib/apt/lists/*

# Install tidyverse
RUN /rocker_scripts/install_tidyverse.sh
# Install Shiny Server
RUN /rocker_scripts/install_shiny_server.sh ${SHINY_VERSION}

# Install R Packages
RUN R -e 'install.packages( c( "abind", "acepack", "anytime", "arm", "askpass", "assertthat", "backports", "base64enc", "BDgraph", "BH", "bindr", "bindrcpp", "bit", "bit64", "blob", "boot", "brio", "broom", "bslib", "bsplus", "cachem", "Cairo", "callr", "car", "carData", "caret", "cellranger", "checkmate", "class", "cli", "clipr", "clock", "cluster", "cmprsk", "coda", "codetools", "colorspace", "combinat", "commonmark", "conflicted", "corpcor", "corrplot", "cowplot", "cpp11", "crayon", "crosstalk", "curl", "CVST", "d3Network", "data.table", "DBI", "dbplyr", "ddalpha", "deldir", "DEoptimR", "desc", "dhglm", "diagram", "diffobj", "digest", "dimRed", "dplyr", "dplyrAssist", "DRR", "DT", "dtplyr", "e1071", "editData", "effects", "ellipsis", "ERSA", "estimability", "evaluate", "exactRankTests", "fansi", "farver", "fastmap", "fdrtool", "fontawesome", "forcats", "foreach", "foreign", "Formula", "frailtyHL", "fs", "future", "future.apply", "gargle", "GauPro", "generics", "geometry", "GGally", "ggcorrplot", "ggforce", "ggfortify", "ggplot2", "ggpubr", "ggraph", "ggrepel", "ggsci", "ggsignif", "ggtext", "glasso", "globals", "glue", "gmp", "googledrive", "googlesheets4", "gower", "graphlayouts", "gridExtra", "gridtext", "gtable", "gtools", "hardhat", "haven", "hexbin", "hglm", "hglm.data", "highr", "Hmisc", "hms", "htmlTable", "htmltools", "htmlwidgets", "httpuv", "httr", "huge", "ids", "igraph", "insight", "interp", "ipred", "isoband", "iterators", "jpeg", "jquerylib", "jsonlite", "kernlab", "KernSmooth", "km.ci", "KMsurv", "knitr", "kutils", "labeling", "later", "lattice", "latticeExtra", "lava", "lavaan", "lazyeval", "lbfgs", "leaps", "lifecycle", "linprog", "lisrelToR", "listenv", "lme4", "lmtest", "lpSolve", "lubridate", "magic", "magrittr", "maptools", "margins", "markdown", "MASS", "Matrix", "matrixcalc", "MatrixModels", "maxstat", "mdhglm", "memoise", "mgcv", "mi", "mime", "miniUI", "minqa", "mitools", "mnormt", "ModelMetrics", "modelr", "munsell", "mvtnorm", "nlme", "nloptr", "nnet", "nortest", "numDeriv", "OpenMx", "openssl", "openxlsx", "parallelly", "pbapply", "pbivnorm", "pbkrtest", "pillar", "pkgconfig", "pkgload", "plogr", "plotly", "pls", "plyr", "png", "polyclip", "polynom", "praise", "prediction", "prettyunits", "pROC", "processx", "prodlim", "progress", "progressr", "promises", "proxy", "ps", "psych", "purrr", "qgraph", "quadprog", "QuantPsyc", "quantreg", "R6", "ragg", "rappdirs", "Rcmdr", "RcmdrMisc", "RColorBrewer", "Rcpp", "RcppArmadillo", "RcppEigen", "RcppParallel", "RcppProgress", "RcppRoll", "readr", "readstata13", "readxl", "recipes", "registry", "regsem", "relimp", "rematch", "rematch2", "reprex", "reshape", "reshape2", "rintrojs", "rio", "rJava", "rjson", "rlang", "rmarkdown", "robustbase", "rockchalk", "ROI", "rpart", "rpf", "rprojroot", "Rsolnp", "rstatix", "rstudioapi", "RUnit", "rvest", "sandwich", "sass", "scales", "selectr", "sem", "semPlot", "semTools", "sfsmisc", "shape", "shiny", "shinyAce", "shinyBS", "shinycssloaders", "shinydashboard", "shinythemes", "shinyWidgets", "slam", "sourcetools", "sp", "spaMM", "SparseM", "spatial", "SQUAREM", "StanHeaders", "stringi", "stringr", "survey", "survival", "survminer", "survMisc", "sys", "systemfonts", "tcltk2", "testthat", "textshaping", "tibble", "tictoc", "tidygraph", "tidyr", "tidyselect", "tidyverse", "timechange", "timeDate", "tinytex", "truncnorm", "tweenr", "tzdb", "utf8", "uuid", "vctrs", "viridis", "viridisLite", "vroom", "waldo", "whisker", "withr", "xfun", "xlsx", "xlsxjars", "XML", "xml2", "xtable", "yaml", "zip", "zoo" ), repos=c( "https://cran.yu.ac.kr/", "https://healthstat.snu.ac.kr/CRAN/", "https://cran.r-project.org/","https://cran.rstudio.com/" ), Ncpus = 16)'

# Copy Source
RUN rm -rf /srv/shiny-server/*
COPY srv/shiny-server /srv/shiny-server/

WORKDIR /srv/shiny-server/

EXPOSE 8787 3838
# run app
USER shiny
CMD ["/usr/bin/shiny-server"]

