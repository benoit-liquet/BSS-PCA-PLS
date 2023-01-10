# Best Subset Selection for PCA and PLS models

Here we provide a comprehensible detailed vignette to run Best Subset Selection for PCA and PLS models

We illustrate the usage of our approach on two real dataset: _multidrug_ and _eQTL\_Hopx_.

This vignette reproduces details result of our application section **Case Studies** from our submitted article:

_Best Subset Selection for Principal Components Analysis and Partial Least Square models using Continuous Optimization_ B. Liquet, S. Moka and S. Muller. (2023).



## Getting Started

In this vignette we used the following R packages

```
library(mixOmics) # this package is available from Bioconductor
llibrary(rlist)
library(R.utils)
library(fields)
library(xtable)
```

## Section 3.2 BSS for PLS model with univariate response

We provide [here](https://github.com/benoit-liquet/BSS-PCA-PLS/blob/main/Section3/Loss_landscape_PLS1.md) the code to illustrate the workings of our continuous optimization method for an example data with $p = 2$.



## Section 5 code for simualtion study [here](https://github.com/benoit-liquet/BSS-PCA-PLS/blob/main/Section-5/)

   - section 5.2: effect of noise [here](https://github.com/benoit-liquet/BSS-PCA-PLS/blob/main/Section-5/Comparison_sPLS_cPLS_TRUE_noise.R)
   - section 5.3: effect of sample size [here](https://github.com/benoit-liquet/BSS-PCA-PLS/blob/main/Section-5/Comparison_sPLS_cPLS_TRUE_sample_size.R)
   - section 5.4: effect of the sparsity [here](https://github.com/benoit-liquet/BSS-PCA-PLS/blob/main/Section-5/Comparison_sPLS_cPLS_TRUE_sparsity.R)
   - section 5.5: effect of dimension $p$ [here](https://github.com/benoit-liquet/BSS-PCA-PLS/blob/main/Section-5/Comparison_sPLS_cPLS_high_p.R)

## Section 6.1 Illustration of Best Subset Selection for PCA

- This analysis is presented [here](https://github.com/benoit-liquet/BSS-PCA-PLS/blob/main/Section-6.1/Vignette_PCA_BSS.md)
 

## Section 6.2 Illustration of Best Subset Selection for PLS2 

- This analysis is presented [here](https://github.com/benoit-liquet/BSS-PCA-PLS/blob/main/Section-6.2/Vignette_PLS2_BSS.md)

