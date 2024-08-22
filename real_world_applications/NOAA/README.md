# NOAA precipitation data application

To replicate the results of the NOAA precipitation data application proceed as follows.
The results of steps 2 and 3 are saved as .rds-files, which are then used in step 4 to generate the tables and figures.

1. Run the file `root\data\NOAA\download_NOAA_data.R` to download and store the NOAA data set.
2. Conduct parameter estimation and prediction with the different methods by running the following scripts:
    - `CholeskyVL.R`: Vecchia-Laplace approximations with Cholesky-based calculations
    - `IterativeVL.R`: Vecchia-Laplace approximations with iterative methods
    - `FITC.R`: Modified predictive process approximation
    - `GPVecchia.R`: Approach of Zilber and Katzfuss [2021]
3. Run `eval_runtime_nll.R` to measure the runtime for computing the marginal log-likelihood with the different methods.
4. Run `show_results.R` to generate the prediction maps and result tables from Appendix A.12.