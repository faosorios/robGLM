# A robust approach for generalized linear models 

[![R](https://img.shields.io/badge/Made%20with-R%20under%20development-success)](https://cran.r-project.org/)
[![fastmatrix](https://img.shields.io/badge/Using%20robustbase-0.95--0-important)](https://cran.r-project.org/package=robustbase)

Supplementary material to **A robust approach for generalized linear models based on maximum Lq-likelihood procedure** by Felipe Osorio, Manuel Galea and Patricia Gimenez

Code written by: Felipe Osorio

Correspondence author: Felipe Osorio, Email: faosorios.stat@gmail.com

Code tested on:
- R under development (2018-02-21 r74285), running Linux Mint 18.3 (64 bits)
- R version 3.3.0, running OS X 10.13.4 (64 bits)

Attached base packages: stats, graphics, grDevices, utils, datasets, methods, base

Other attached packages (updated): robustbase-0.95

CONTENTS:
- case_study/case_Finney.R: R commands for the analysis of Finney's dataset.
- case_study/case_leukemia.R: R commands for the analysis of Leukemia dataset.
- case_study/suppl_Finney.R: R commands for additional result reported in the Supplementary Material
- code/envelope.R: R functions for computation of QQ-plots with simulated envelopes.
- code/qopt.R: R function for the selection of distortion parameter 'q'.
- code/robLogistic.R: R functions for robust estimation in the logistic regression model (used to analyse Finney's data).
- code/robLogLin.R: R functions for robust estimation in the log-linear regression model (used in simulation experiment).
- data/finney.csv: transient vasoconstriction (or Finney's) dataset in CSV format.
- data/finney.rda: transient vasoconstriction (or Finney's) dataset in RDA format.
- data/leukemia.rda: Leukemia dataset in RDA format.
- simul/simul_Poi.R: R function to setting the model, data contamination and fitting log-linear models.
- simul/simulation.R: R commands to perform the Monte Carlo simulation study.
- README.md: this file.
