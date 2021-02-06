## Simulation testing survey geostatistical index standardization

### Required R packages to install:

- SimSurvey: https://github.com/PaulRegular/SimSurvey
- sdmTMB: https://github.com/pbs-assess/sdmTMB

### General examples:

- [`example.Rmd`](example.Rmd): a complete walk through of a basic example
- [`example-rep.R`](example-rep.R): a basic example of simulation testing with replicates

### Examples for specific research topics:

- [`covariate-index.R`](covariate-index.R): demonstrates the effect of including a depth covariates when calculating a model-based index

- [`stitching-index.R`](stitching-index.R): demonstrates the effect of attempting to "stitch" contiguous surveys with different catchabilities together

- [`coverage-index.R`](coverage-index.R): demonstrates the effect of missing spatial survey coverage in some years
