# Simulation studies

Each script performs its own simulation study and at the end of the script a plot function is included to display the results. Simulated data is generated using the function `make_data()` in the file `root\data\simulated\make_data.R`.

* `bias_analysis.R`

  Covariance parameter estimation for different sample sizes using iterative methods for Vecchia-Laplace approximations as shown in Figure 1, Section 1.

* `compare_preconditioner_nll.R`

  Preconditioner comparison (VADU and LRAC) with regard to runtime and variance of log-marginal likelihoods for different numbers of random vectors.
  Set the range (`rho_true`) to 0.05 (Figure 1), 0.01 (Figure 2), and 0.25 (Figure 3) to replicate the results in Appendix A.8.

* `estimation_prediction_n2e4.R`

  Covariance parameter estimation and prediction on 100 simulated data sets with sample size n=n_p=20'000 using Iterative-VL, Cholesky-VL, and GPVecchia.
    - Evaluation of the relative differences in the log-likelihood compared to Cholesky-based calculations as shown in Figure 6, Section 4.6.
    - Evaluation of covariance parameter estimates as shown in Figure 4, Section 4.4, and Table 1, Appendix A.10.
    - Evaluation of the prediction accuracy as shown in Figure 5, Section 4.5, and Table 2, Appendix A.10.
    - Evaluation of the runtime for estimation as shown in Table 1, Section 4.7.

* `eval_predictive_variances.R`

  Compare the accuracy and runtime for simulation-based predictive variances relative to Cholesky-based computations as shown in Figure 5, Appendix A.9.
  
* `eval_runtime_nll.R`

  Runtime comparison to evaluate the marginal log-likelihood for different sample sizes as shown in Figure 7, Section 4.7.
  
* `NN_nll_estimation_prediction_n5e3.R`

  Varying number of nearest neighbors for the Vecchia-Laplace approximation and a sample size of n=n_p=5000. 
  Evaluation of the marginal likelihood, prediction accuracy, and runtime for the marginal log-likelihood for Iterative-VL and GPVecchia relative to Laplace as shown in the upper part of Figure 8, Section 4.8.
  
* `NN_nll_prediction_n5e4.R`

  Varying number of nearest neighbors for the Vecchia-Laplace approximation and a sample size of  n=n_p=20'000.
  Evaluation of the marginal likelihood, prediction accuracy, and runtime for the marginal log-likelihood for Iterative-VL and GPVecchia as shown in the lower part of Figure 8, Section 4.8.
  
* `rank_LRAC.R`

  Compare different ranks for the LRAC preconditioner with regard to runtime and variance of log-marginal likelihoods as shown in Figure 4, Appendix A.8.
  