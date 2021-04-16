# KKT_1D_GTV

This repo implements the algorithm in the paper [A unified approach for a 1D generalized total variation problem. Mathematical Programming, February, 2021.](http://link.springer.com/article/10.1007/s10107-021-01633-2)

The problem to be solved is of the form:

\min_{x_1,\ldots, x_n} \sum_{i=1}^n f_i(x_i) + \sum_{i=1}^{n-1} h_i(x_i - x_{i+1}).

Both f_i(x_i) and h_i(x_i - x_{i+1}) functions are arbitrary convex functions.

## Reference

Please cite the paper if you use the algorithm.

**Lu, C., Hochbaum, D.S.** A unified approach for a 1D generalized total variation problem. *Math. Program.* (2021). https://doi.org/10.1007/s10107-021-01633-2
