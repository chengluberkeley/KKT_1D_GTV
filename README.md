# KKT_1D_GTV

This repo implements the algorithm in the paper [A unified approach for a 1D generalized total variation problem. Mathematical Programming, February, 2021.](http://link.springer.com/article/10.1007/s10107-021-01633-2)

The problem to be solved is of the form:

<img src="https://render.githubusercontent.com/render/math?math=\text{(1D-GTV)}\quad\min_{x_1,\ldots, x_n}\sum_{i=1}^n f_i(x_i) %2B \sum_{i=1}^{n-1} h_i(x_i - x_{i %2B 1})">

Both <img src="https://render.githubusercontent.com/render/math?math=f_i(x_i)"> and <img src="https://render.githubusercontent.com/render/math?math=h_i(x_i - x_{i %2B 1})"> functions are arbitrary convex functions.

## Special cases
Below are some special cases of 1D-GTV that are implemented in the code:
#### L2-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n}\frac{1}{2}\sum_{i=1}^nc_i(x_i - a_i)^2 %2B \sum_{i=1}^{n-1}c_{i,i%2B 1}|x_i - x_{i %2B 1}|">

#### piecewise-quadratic-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \sum_{i=1}^nf^{pq}_i(x_i) %2B \sum_{i=1}^{n-1}c_{i,i%2B1}|x_i - x_{i %2B 1}|">

Each <img src="https://render.githubusercontent.com/render/math?math=f^{pq}_i(x_i)"> is a convex piecewise quadratic function.

#### L2-Tikhonov
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n}\frac{1}{2}\sum_{i=1}^nc_i(x_i - a_i)^2 %2B \frac{1}{2}\sum_{i=1}^{n-1}c_{i,i%2B 1}(x_i - x_{i %2B 1})^2">

#### L1-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n}\sum_{i=1}^nc_i|x_i - a_i| %2B \sum_{i=1}^{n-1}c_{i,i%2B 1}|x_i - x_{i %2B 1}|">

#### piecewise-linear-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \sum_{i=1}^nf^{pl}_i(x_i) %2B \sum_{i=1}^{n-1}c_{i,i%2B1}|x_i - x_{i %2B 1}|">

Each <img src="https://render.githubusercontent.com/render/math?math=f^{pl}_i(x_i)"> is a convex piecewise linear function.

#### Lp-Lq
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n}\frac{1}{p}\sum_{i=1}^nc_i|x_i - a_i|^p %2B \frac{1}{q}\sum_{i=1}^{n-1}c_{i,i%2B1}|x_i - x_{i %2B 1}|^q">

#### Huber-TV
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \sum_{i=1}^nc_i\rho_{k_i}(x_i-a_i) %2B \sum_{i=1}^{n-1}c_{i,i%2B1}|x_i - x_{i %2B 1}|">

Each <img src="https://render.githubusercontent.com/render/math?math=\rho_{k_i}(x_i-a_i)"> is a [Huber](https://en.wikipedia.org/wiki/Huber_loss) loss function defined as:
<img src="https://render.githubusercontent.com/render/math?math=\rho_{k_i}(x_i-a_i) = \begin{cases}\frac{1}{2}(x_i-a_i)^2,\quad |x_i-a_i| \leq k_i \\
k_i|x_i - a_i| - \frac{1}{2}k^2_i,\quad \text{otherwise}
\end{cases}">

#### L2-Huber
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \frac{1}{2}\sum_{i=1}^nc_i(x_i - a_i)^2 %2B \sum_{i=1}^{n-1}c_{i,i%2B1}\rho_{k_{i,i%2B1}}(x_i - x_{i %2B 1})">

#### Linear-L2 (Laplacian on a path)
<img src="https://render.githubusercontent.com/render/math?math=\min_{x_1,\ldots,x_n} \sum_{i=1}^nc_ix_i %2B \frac{1}{2}\sum_{i=1}^{n-1}(x_i - x_{i %2B 1})^2">

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
cmake .. && make . -j5
```

## Solve your own problem
To solve (1D-GTV) problem of your own <img src="https://render.githubusercontent.com/render/math?math=f_i(x_i)"> and <img src="https://render.githubusercontent.com/render/math?math=h_i(x_i - x_{i %2B 1})"> functions, simply update the functions `compDrvt(...)` and `compSepInv(...)` in [kkt.hpp](KKT/kkt.hpp) to compute the derivatives and the inverses of derivatives for your <img src="https://render.githubusercontent.com/render/math?math=f_i(x_i)"> and <img src="https://render.githubusercontent.com/render/math?math=h_i(x_i - x_{i %2B 1})"> functions respectively.

## Reference

Please cite the paper if you use the algorithm.

**Lu, C., Hochbaum, D.S.** A unified approach for a 1D generalized total variation problem. *Math. Program.* (2021). https://doi.org/10.1007/s10107-021-01633-2
