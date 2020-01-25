runAll compares different implementations of our spectral bundle method for minimizing 
\lambda_max problem. The function is defined at maxEig.m. 
See SB_ls, SBQ_NKKT and  SBQ_CVX for further details. An example is

runAll(1, 1000, 0.01,40,50)

with arguments:

genprob = 1; generates a random max_eig problme with matrix dim = 50 and # of constraints=40;
maxit = 1000;
eigTol = 0.01, is being used to guess the optimal multiplicity of the lam_max.

