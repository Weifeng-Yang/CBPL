# CBPL
This package contains code for the Controlled Block Proximal Linearized Method (CBPL) proposed in our paper[<sup>1</sup>](#refer-id) for solving the graph-regularized sparse non-negative matrix factorization with $\ell_0$-constraints ($\ell_0$-GSNMF) problem, the graph-regularized sparse non-negative CP decomposition with $\ell_0$-constraints ($\ell_0$-GSNCP) problem, and the graph-regularized sparse non-negative independent-core matrix tri-factorization with $\ell_0$-constraints ($\ell_0$-GSNMTF) problem.

## Matlab code
"L0GSNMF" is the $\ell_0$-GSNMF code, "L0GSNCP" is the $\ell_0$-GSNCP code, and "L0GSNMTF" is the $\ell_0$-GSNMTF code.

## How to run our code
For "L0GSNMF", just run '[main_Run_me.m](L0GSNMF/main_Run_me.m)'.

For "L0GSNCP", just run '[main_Run_me.m](L0GSNCP/main_Run_me.m)'.

For "L0GSNMTF", just run '[main_Run_me.m](L0GSNMTF/main_Run_me.m)'.

## Data
This code has built-in the data mentioned in our paper[<sup>1</sup>](#refer-id).

In each problem folder, "Data" contains the input datasets, and "Data_run" contains the saved experimental results, including objective function values, running times, residuals and factorization results.

## Reference
<div id="refer-id"></div>
[1] Controlled Block Proximal Linearized Method for Nonconvex and Nonsmooth Optimization. Weifeng Yang. 

