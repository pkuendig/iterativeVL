# MODIS column water vapor data application

To replicate the results of the MODIS column water vapor data application proceed as follows.
The results of steps 3 and 4 are saved as .rds-files, which are then used in step 5 to generate the tables and figures.

1. Run `root\data\MYD05\read_hdf.R` to read the data in the .hdf-file.
2. Run `root\data\MYD05\sample_MYD05.R` to sample the training and test data sets.
3. Conduct parameter estimation and prediction with the different methods by running the following scripts:
    - `CholeskyVL.R`: Vecchia-Laplace approximations with Cholesky-based calculations
    - `IterativeVL.R`: Vecchia-Laplace approximations with iterative methods
    - `FITC.R`: Modified predictive process approximation
    - `GPVecchia.R`: Approach of Zilber and Katzfuss [2021]
    - `IterativeVL_n25e4.R` (optional): Estimation on the large training data set for Vecchia-Laplace approximations with iterative methods (Table 4 in Appendix A.11).
4. Run `eval_runtime_nll.R` to measure the runtime for computing the marginal log-likelihood with the different methods.
5. Run `show_results.R` to generate the prediction maps and result tables from Section 5 and Appendix A.11.