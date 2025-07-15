# Contents
An overall description of the files in the code folder, which we used to run our analysis is found here.
# System requirements
We run the analysis in the following R environment and package versions:
# Installation guide
To install R, please follow the instructions found here.

To install the R packages needed for our analysis you can run the following:

install.packages(
    c(
      "readr", "dplyr", "tidyr", "purrr", "lubridate",
      "uwot", "igraph", "mvtnorm", "rio", "stringr", 
      "ggh4x", "patchwork", "ggraph", "mclust", 
      "archetypes", "dbscan", "survival", "ggflowchart", 
      "reshape2", "scales", "lme4", "lmerTest", "meta", 
      "broom", "broom.mixed", "MGMM", "dcurves", "Rtsne", 
      "kernlab", "ClustOfVar", "GGally", "ggdensity", "glmnet"
    )
)
On a typical desktop computer this should take around 1 hour.

# Reproduction
Due to data access restrictions, we cannot include the data necessary to reproduce all the quantitative results in the manuscript. However, it is possible to run the same analysis in a new dataset provided it has the same format as we used.
