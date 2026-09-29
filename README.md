# Calcium Atlas (Ca²⁺tlas): analysis code

Code for:

Župunski M, Brugman R, Zauser M, Lampou K, Waadt R, Guichard M, Krebs M, Schumacher K, Lamar R, Pahle J, Grossmann G (2026). *A calcium-fingerprint atlas for inferring the signaling modes of complex elicitors, applied to a humic biostimulant.*

## Contents

`R/Crestline_Ca2tlas_analysis.R` reproduces every statistic, table and quantitative plot in the paper from the Supplemental Datasets (sheets "dataset 1" to "dataset 10"):

- descriptor comparisons across the 19 reference elicitors (Fig. S4)
- the Ca²⁺tlas PCA and within-treatment correlations (Fig. 5)
- projection of the biostimulant BC4 into the atlas space (Fig. 6b)
- class-conditional multivariate normal (MVN) classifier: posterior membership under frequency and uniform priors, Mahalanobis distances, typicality and kNN (Fig. 6c, S6; Table S3)
- stratified 5-fold cross-validation (Fig. S6d, Table S2) and leave-one-elicitor-out validation (Table S4)
- biosensor and imaging-platform comparisons (Fig. S5)
- reporter-line validation and BC4 responses of *moca1* and *cerk1-2* (Fig. 7, S8, S9, S11; Table S5)
- primary root growth models (Fig. 7e, S10b)
- intensity traces (Fig. 1, S1, S7, S9c)

The script ends with a reproduction check that compares each number quoted in the manuscript with the value it computes.

## How to run

1. Download the Supplemental Datasets workbook and save it as `Supplemental_Datasets.xlsx`.
2. Open the script in R (tested with R 4.3) and set `setwd()` in Block 0 to the folder that holds the workbook.
3. Run Block 0 and Block 1, then any other block. Each block header states what it needs.

Missing packages are installed automatically in Block 0: readxl, dplyr, tidyr, ggplot2, FactoMineR, factoextra, corpcor, mvtnorm, FSA, multcompView, Hmisc, emmeans, car.

All results are written to the folder `Crestline_outputs`.

## Fiji macros

The Fiji macros that generate kymographs and crestplots from the image stacks: [link to be added].

## License

MIT (see LICENSE).
