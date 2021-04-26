# KKT_1D_GTV

This repo implements the algorithm in the paper [A unified approach for a 1D generalized total variation problem. Mathematical Programming, February, 2021.](http://link.springer.com/article/10.1007/s10107-021-01633-2)

The problem to be solved is of the form:

<img src="https://render.githubusercontent.com/render/math?math=\text{(1D-GTV)}\quad\min_{x_1,\ldots, x_n}\sum_{i=1}^n f_i(x_i) %2B \sum_{i=1}^{n-1} h_i(x_i - x_{i %2B 1})">

Both <img src="https://render.githubusercontent.com/render/math?math=f_i(x_i)"> and <img src="https://render.githubusercontent.com/render/math?math=h_i(x_i - x_{i %2B 1})"> functions are arbitrary convex functions.

## Special cases
Below are some special cases of 1D-GTV that are implemented in the code:
#### L2-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n}\frac{1}{2}\sum_{i=1}^nc_i(x_i - a_i)^2 %2B \sum_{i=1}^{n-1}c_{i,i%2B 1}|x_i - x_{i %2B 1}|">

Sample run on data generated according to [1] (non-weighted):
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 | 10000000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 2 | 8.4 | 60.6 | 2647.4 | 8952.8 |
| Runtime std (ms) | 0 | 0 | 0 | 0.8 | 3.07246 | 93.8266 | 278.203 |

[1] Condat, L.: A direct algorithm for 1-D total variation denoising. IEEE Signal Process. Lett. 20(11), 1054--1057 (2013).

Sample run on data generated according to [2] (weighted):
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 | 10000000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 0 | 3 | 28 | 271.6 | 2718 |
| Runtime std (ms) | 0 | 0 | 0 | 0 | 1.26491 | 11.6722 | 30.9839 |

[2] Barbero, Á., Sra, S.: Modular proximal optimization for multidimensional total-variation regularization. J. Mach. Learn. Res. 19(1), 2232--2313 (2018).

#### piecewise-quadratic-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \sum_{i=1}^nf^{pq}_i(x_i) %2B \sum_{i=1}^{n-1}c_{i,i%2B1}|x_i - x_{i %2B 1}|">

Each <img src="https://render.githubusercontent.com/render/math?math=f^{pq}_i(x_i)"> is a convex piecewise quadratic function.

Sample run times for increasing `n`:
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 |
| --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 2 | 18.6 | 198.4 | 2094 |
| Runtime std (ms) | 0 | 0 | 0 | 1.0198 | 6.97424 | 40.5364 |

Sample run times for increasing number of breakpoints per <img src="https://render.githubusercontent.com/render/math?math=f^{pq}_i(x_i)">:
| #breakpoints | 100 | 200 | 300 | 400 | 500 | 600 |
| --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 198.2 | 241.2 | 279 | 290.8 | 300 | 313.4 |
| Runtime std (ms) | 1.16619 | 11.0164 | 7.79744 | 5.49181 | 8.50882 | 6.74092 |
`n` is fixed to 100000.

#### L2-Tikhonov
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n}\frac{1}{2}\sum_{i=1}^nc_i(x_i - a_i)^2 %2B \frac{1}{2}\sum_{i=1}^{n-1}c_{i,i%2B 1}(x_i - x_{i %2B 1})^2">

Sample run times:
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 | 10000000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 1 | 15.2 | 160.8 | 1550.6 | 15403.2 |
| Runtime std (ms) | 0 | 0 | 0 | 1.16619 | 9.36803 | 39.0056 | 175.991 |

#### L1-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n}\sum_{i=1}^nc_i|x_i - a_i| %2B \sum_{i=1}^{n-1}c_{i,i%2B 1}|x_i - x_{i %2B 1}|">

Sample run times:
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 | 10000000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 0 | 4.2 | 54.2 | 497.2 | 4928.6 |
| Runtime std (ms) | 0 | 0 | 0 | 0.4 | 3.18748 | 9.08625 | 174.691 |

#### piecewise-linear-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \sum_{i=1}^nf^{pl}_i(x_i) %2B \sum_{i=1}^{n-1}c_{i,i%2B1}|x_i - x_{i %2B 1}|">

Each <img src="https://render.githubusercontent.com/render/math?math=f^{pl}_i(x_i)"> is a convex piecewise linear function.

Sample run times for increasing `n`:
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 |
| --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 1 | 17 | 177 | 1772.6 |
| Runtime std (ms) | 0 | 0 | 0 | 0 | 2.44949 | 37.4945 |

Sample run times for increasing number of breakpoints per <img src="https://render.githubusercontent.com/render/math?math=f^{pl}_i(x_i)">:
| #breakpoints | 100 | 200 | 300 | 400 | 500 | 600 |
| --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 171.6 | 226.8 | 255.6 | 271.8 | 282.4 | 296 |
| Runtime std (ms) | 1.35647 | 3.81576 | 1.49666 | 2.92575 | 2.41661 | 1.89737 |
`n` is fixed to 100000.

#### Lp-Lq
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n}\frac{1}{p}\sum_{i=1}^nc_i|x_i - a_i|^p %2B \frac{1}{q}\sum_{i=1}^{n-1}c_{i,i%2B1}|x_i - x_{i %2B 1}|^q">

Sample run times of `p = q = 4`:
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 | 10000000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 2 | 21.4 | 213.4 | 1997 | 20344.4 |
| Runtime std (ms) | 0 | 0 | 0 | 0.489898 | 18.4781 | 18.868 | 396.665 |

#### Huber-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \sum_{i=1}^nc_i\rho_{k_i}(x_i-a_i) %2B \sum_{i=1}^{n-1}c_{i,i%2B1}|x_i - x_{i %2B 1}|">

Each <img src="https://render.githubusercontent.com/render/math?math=\rho_{k_i}(x_i-a_i)"> is a [Huber](https://en.wikipedia.org/wiki/Huber_loss) loss function defined as:
<img src="https://render.githubusercontent.com/render/math?math=\rho_{k_i}(x_i-a_i) = \begin{cases}\frac{1}{2}(x_i-a_i)^2,\quad |x_i-a_i| \leq k_i \\
k_i|x_i - a_i| - \frac{1}{2}k^2_i,\quad \text{otherwise}
\end{cases}">

Sample run times:
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 | 10000000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 0 | 3.6 | 45.8 | 435.2 | 4077.8 |
| Runtime std (ms) | 0 | 0 | 0 | 0.489898 | 1.16619 | 27.0215 | 60.562 |

#### L2-Huber
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \frac{1}{2}\sum_{i=1}^nc_i(x_i - a_i)^2 %2B \sum_{i=1}^{n-1}c_{i,i%2B1}\rho_{k_{i,i%2B1}}(x_i - x_{i %2B 1})">

Sample run times:
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 | 10000000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 0 | 8.4 | 81.4 | 826.8 | 8143.4 |
| Runtime std (ms) | 0 | 0 | 0 | 0.8 | 3.38231 | 10.1863 | 137.283 |

#### Linear-L2 (Laplacian on a path)
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \sum_{i=1}^nc_ix_i %2B \frac{1}{2}\sum_{i=1}^{n-1}(x_i - x_{i %2B 1})^2">

Sample run times:
| Problem scale (n) | 10 | 100 | 1000 | 10000 | 100000 | 1000000 | 10000000 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Average runtime (ms) | 0 | 0 | 0 | 0 | 0 | 1.4 | 17.2 |
| Runtime std (ms) | 0 | 0 | 0 | 0 | 0 | 0.8 | 9.92774 |

For all the above experimental results:
* Each average is taken by 5 runs of randomly (if not specified) generated data and coefficients.
* All runs are on a MacBook (Big Sur) of 2.4 GHz 8-Core Intel Core i9, 32 GB 2667 MHz DDR4 memory.

## Run prebuilt executable
For a sample run of all special case problems implemented:
```
cd bin
./run.sh
```
It will output runtime profiles in the `/tmp/` folder.

To see the full runtime arguments,
```
cd bin
./kkt_main
```

## Build from the source
This project is a `cmake` project. To build from the source:
```
mkdir build && cd build
cmake .. && make -j5
```

## Solve your own problem
To solve (1D-GTV) problem of your own <img src="https://render.githubusercontent.com/render/math?math=f_i(x_i)"> and <img src="https://render.githubusercontent.com/render/math?math=h_i(x_i - x_{i %2B 1})"> functions, simply update the functions `compDrvt(...)` and `compSepInv(...)` in [kkt.hpp](KKT/kkt.hpp) to compute the derivatives and the inverses of derivatives for your <img src="https://render.githubusercontent.com/render/math?math=f_i(x_i)"> and <img src="https://render.githubusercontent.com/render/math?math=h_i(x_i - x_{i %2B 1})"> functions respectively.

## Reference

Please cite the paper if you use the algorithm.

**Lu, C., Hochbaum, D.S.** A unified approach for a 1D generalized total variation problem. *Math. Program.* (2021). https://doi.org/10.1007/s10107-021-01633-2
