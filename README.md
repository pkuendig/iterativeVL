# Iterative Methods for Vecchia-Laplace Approximations
This repository contains the R code for the simulations in the paper "Iterative Methods for Vecchia-Laplace Approximations for Latent Gaussian Process Models".

The iterative methods for the Vecchia-Laplace approximation are implemented in the package GPBoost and available at https://github.com/fabsig/GPBoost. Release v1.2.5 is recommended for reproducibility of the results in the paper.

Further comments:

- The provided MYD05 data set (MYD05_L2.A2019087.1345.061.2019088165734.hdf) is from the NASA Earthdata website: https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MYD05_L2#overview. This website provides more information about the data and a link to the data archive.
- To replicate the real-world application, first process the data using the provided R files in the MYD05 folder.
- We do not provide the simulations in this repository that use Lanczos-based prediction or preconditioners other than the VADU and LRAC, since these extensions are not included in the public version of GPBoost.
