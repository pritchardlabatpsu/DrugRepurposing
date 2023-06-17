# Drug Repurposing

This project aims to predict clinically relevant in vitro drug concentrations from pharmacokinetic parameters and drug-protein interaction assays.

The file structure is as follows,
```
.
├── Figures            # code to generate .ai figure files from raw data
│   ├── Figure1        # meta-analysis of clinical trial success across TKIs and preclinical analysis of imatinib and erlotinib
│   ├── Figure2        # heatmaps of various metrics of in vitro dose response curves (IC50-fold change, free alpha, effective alpha, etc) and associated ROC curves
│   ├── Figure3        # mapping of closed-form solution of analytical drug binding model to parameter space and analysis of protein titration experiments
│   ├── Figure4        # logistic models and their evaluation, including the final effective C_ave model (LogisticRegressionModels.rds)
│   ├── Figure5        # results of logistic effective C_ave model for EGFR and KIT inhibitors
│   ├── Figure6        # synthesis of clinical and preclinical meta-analysis for imatinib
├── ReferenceFiles     # where raw data is stored
├── Supplement         # code and resources to generate manuscript supplement
│   ├── SuppAppendices # code to map solutions of analytical off-targeting drug binding model
│   ├── SuppFigures    # collection of dendrogram, GR50, serum protein titration, and free C_ave analyses
│   ├── SuppMethods    # description of clinical and preclinical meta-analysis methods
│   ├── SuppTables     # collection of supplemental tables used in analysis
├── WebApp             # raw code for R shiny app for custom serum-shift experiments

```

## Installation
After cloning the repository, install the following required R packages

```
# Package names
packages = c("rstudioapi", "ggplot2", "RColorBrewer", "scales", "reshape2", "pROC", "viridis", "dr4pl")

# Install packages not yet installed
installed_packages = packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}
```

## Data

The repo here contains all data and code necessary to generate the logistic classifiers and analyze the results. Raw data can be found in the ReferenceFiles folder.

## Manuscript

The code used to generate figures in the manuscript are organized in the Figures folder.