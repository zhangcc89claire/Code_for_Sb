# Unraveling NEK4 as a Potential Drug Target in Schizophrenia and Bipolar I Disorder: Analysis Pipeline

This repository contains the analysis code for the paper "Unraveling NEK4 as a Potential Drug Target in Schizophrenia and Bipolar I Disorder: A Proteomic and Genomic Approach". The comprehensive analysis pipeline encompasses methodologies from Proteome-Wide Association Studies (PWAS), Transcriptome-Wide Association Studies (TWAS), Mendelian Randomization (MR) and CELLEX analysis, providing a robust approach to understanding the genetic underpinnings and potential therapeutic targets for Schizophrenia and Bipolar I Disorder.

## Prerequisites

Before starting with the analysis, ensure that the following software and libraries are installed and configured on your system:

- R and the following R packages: `FUSION`, `TwoSampleMR`
- Python and the following Python packages: `numpy`, `pandas`, `cellex`

## Analysis Steps

### Step 1: PWAS

For the Proteome-Wide Association Study (PWAS), execute the `FUSION.assoc_test.R` script using prepared data for all the chromosomes. Below, we provide an example using chromosome 22:
  
```bash
Rscript FUSION.assoc_test.R \
  --sumstats gwas.sumstats \
  --weights ./WEIGHTS/gene.pos \
  --weights_dir ./WEIGHTS/ \
  --ref_ld_chr ./LDREF/1000G.EUR. \
  --chr 22 \
  --out PGC2.SCZ.chr22.dat



### Step 2: TWAS

For the Transcriptome-Wide Association Study (TWAS), use the same script with RNA position weights. Below, we provide an example using chromosome 22:

Rscript FUSION.assoc_test.R \
--sumstats gwas.sumstats \
--weights ./WEIGHTS/gene.pos \
--weights_dir ./WEIGHTS/ \
--ref_ld_chr ./LDREF/1000G.EUR. \
--chr 22 \
--out PGC2.SCZ.chr.dat


### Step 3: MR
Finally, for the Mendelian Randomization analysis, use the TwoSampleMR package in R to explore the causal relationship.

exposure_dat <- read_excel("genetic_instrument_protein.xlsx")
outcome_dat <- read_xlsx("outcome_SCZ or BDI.xlsx")
dat <- harmonise_data(exposure_dat, outcome_dat)
res <- mr(dat)

### Step 4: CELLEX

To perform the CELLEX analysis for cell-type specific expression insights, use the following Python code.

import numpy as np
import pandas as pd
import cellex



#### Loading data
data = pd.read_csv("./data.csv", index_col=0)
metadata = pd.read_csv("./metadata.csv", index_col=0)

# Initializing and computing ESObject
eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
eso.compute(verbose=True)

# Saving results
eso.results["esmu"].to_csv("mydataset.esmu.csv.gz")
