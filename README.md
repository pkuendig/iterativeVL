# Iterative Methods for Vecchia-Laplace Approximations

This repository contains the R code for the simulation studies and real-world applications in the paper "Iterative Methods for Vecchia-Laplace Approximations for Latent Gaussian Process Models". Iterative methods for Vecchia-Laplace approximations are implemented in the package GPBoost: <https://github.com/fabsig/GPBoost>.

* The `simulation_studies` folder contains the scripts to run the simulations. See `simulation_studies\README.md` for a detailed description of the different simulations.
* The `real_world_applications` folder contains the scripts to run the real-world applications for the MODIS column water vapor data (`real_world_applications\MYD05`) and the NOAA precipitation data (`real_world_applications\NOAA`). The `README.md` file in each subdirectory gives step-by-step instructions on how to run the real-world application.
* The `data` folder contains all the files needed to generate the necessary data.

### Structure of the repository
```
root
│   README.md
│
└───data
│   │
│   └───MYD05
│   │   │   MYD05_L2.A2019087.1345.061.2019088165734.hdf
│   │   │   read_hdf.R
│   │   │   sample_MYD05.R
│   │
│   └───NOAA
│   │   │   download_NOAA_data.R
│   │
│   └───simulated
│       │   make_data.R
│
└───real_world_applications
│   │
│   └───MYD05
│   │   │   README.md
│   │   │   CholeskyVL.R
│   │   │   eval_runtime_nll.R
│   │   │   FITC.R
│   │   │   GPVecchia.R
│   │   │   IterativeVL.R
│   │   │   IterativeVL_n25e4.R
│   │   │   show_results.R 
│   │
│   └───NOAA
│       │   README.md
│       │   CholeskyVL.R
│       │   eval_runtime_nll.R
│       │   FITC.R
│       │   GPVecchia.R
│       │   IterativeVL.R
│       │   show_results.R  
│
└───simulation_studies
    │   README.md
    │   bias_analysis.R
    │   compare_preconditioner_nll.R
    │   estimation_prediction_n2e4.R
    │   eval_predictive_variances.R
    │   eval_runtime_nll.R
    │   NN_nll_estimation_prediction_n5e3.R
    │   NN_nll_prediction_n5e4.R
    │   rank_LRAC.R
```
The provided data set `MYD05_L2.A2019087.1345.061.2019088165734.hdf` is from the NASA Earthdata website (<https://ladsweb.modaps.eosdis.nasa.gov/missions-and-measurements/products/MYD05_L2#overview>). 