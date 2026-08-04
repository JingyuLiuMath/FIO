# FIO

## Overview

This repository contains MATLAB research code for discrete Fourier integral operators (FIOs).

## Discrete Fourier Integral Operators

The $d$-dimensional discrete FIO takes the form

$$
u(x) = \sum_{\xi \in \Xi} a(x, \xi)
\exp(2 \pi \imath \phi(x, \xi)) f(\xi),
\qquad x \in X,
$$

where

$$
X = \bigl\{
x_{i} = (i_{1} / n, \dotsc, i_{d} / n)
: 0 \leq i_{1}, \dotsc, i_{d} < n
\bigr\}
$$

and

$$
\Xi = \bigl\{
\xi_{j} = (j_{1}, \dotsc, j_{d})
: -n / 2 \leq j_{1}, \dotsc, j_{d} < n / 2
\bigr\}
$$

are the spatial and frequency grids, respectively. Here, $n$ is a positive even integer.

The amplitude $a(x, \xi)$ is assumed to be smooth in $(x, \xi)$. The phase $\phi(x, \xi)$ is smooth for $\xi \neq 0$ and homogeneous of degree one in $\xi$:

$$
\phi(x, \lambda \xi) = \lambda \phi(x, \xi),
\qquad \lambda > 0.
$$

The number of degrees of freedom is $N = n^d$. In matrix form, the discrete FIO is written as

$$
u = Kf.
$$

## Butterfly Factorization (BF)

### Sparse Matrix Product (SMP) Form

A butterfly factorization approximates $K$ by

$$
K \approx \widetilde{K}
= U^{[L]} G^{[L]} \dotsb G^{[h + 1]} M^{[h]}
H^{[h + 1]} \dotsb H^{[L]} V^{[L]},
$$

where $L = \mathcal{O}(\log n)$ and $h \approx L / 2$. For fixed butterfly ranks, each sparse factor has $\mathcal{O}(N)$ nonzero entries.
Consequently, applying $\widetilde{K}$ or $\widetilde{K}^{*}$ costs $\mathcal{O}(N \log N)$ for fixed dimension. For the FIOs considered here, the factorization can be constructed by interpolation in
$\mathcal{O}(N \log N)$ operations when the interpolation rank is treated as a constant.

This form of BF is implemented by [FastBF](https://github.com/YingzhouLi/FastBF.m), which needs to be installed under the `extern` directory.

Related papers:

- Yingzhou Li, Haizhao Yang, Eileen R. Martin, Kenneth L. Ho, and Lexing Ying, *Butterfly Factorization*, Multiscale Modeling & Simulation, 13 (2015), pp. 714--732.
- Yingzhou Li and Haizhao Yang, *Interpolative butterfly factorization*, SIAM Journal on Scientific Computing, 39 (2017), pp. A503--A531.

### Block Hadamard Product (BHP) Form

Alternatively, the BF can be written as

$$
K
\approx \widetilde{K}
= U^{[L]} \odot G^{[L]} \odot \dotsb \odot G^{[h + 1]}
\odot M^{[h]} \odot H^{[h + 1]} \odot \dotsb \odot H^{[L]}
\odot V^{[L]},
$$

where each factor is an array of small dense blocks. Here, $\odot$ denotes a generalized block Hadamard product: corresponding dense blocks are multiplied as matrices, rather than entrywise, and are
redistributed according to the complementary trees between successive levels. For fixed butterfly ranks and dimension, this representation retains the $\mathcal{O}(N\log N)$ storage and application complexities of BF while avoiding the expanded intermediate vectors arising from global sparse matrix products. It is implemented by the `Butterfly` and `Butterfly2D` classes in this repository.

## Hierarchically Semiseparable (HSS) Matrices

Let $\mathsf{T}$ be a hierarchical partition of an index set $\mathcal{J}$. For the Hermitian matrices considered here, the two defining properties of the HSS representation are low-rank off-diagonal blocks and nested bases. For example, if a node $\tau$ has two children $\alpha_{1}$ and $\alpha_{2}$, then

$$
D_{\tau}
=
\begin{bmatrix}
D_{\alpha_{1}}
& U_{\alpha_{1}}^{\mathrm{big}} B_{\alpha_{2}, \alpha_{1}}^{*}
  U_{\alpha_{2}}^{\mathrm{big}, *} \\
U_{\alpha_{2}}^{\mathrm{big}} B_{\alpha_{2}, \alpha_{1}}
  U_{\alpha_{1}}^{\mathrm{big}, *}
& D_{\alpha_{2}}
\end{bmatrix},
$$

where $D_{\tau} = H(\mathcal{J}_{\tau}, \mathcal{J}_{\tau})$, and

$$
U_{\tau}^{\mathrm{big}}
=
\begin{bmatrix}
U_{\alpha_{1}}^{\mathrm{big}} & 0 \\
0 & U_{\alpha_{2}}^{\mathrm{big}}
\end{bmatrix}
U_{\tau}.
$$

For 1D problems, a binary tree is used. For 2D problems, a quadtree is used.

## Approximate Inversion via HSS Matrices

Assume that $K$ is square and nonsingular. Define the Hermitian positive definite normal matrix

$$
G = K^{*} K.
$$

The inverse of $K$ is given by

$$
K^{-1} = G^{-1} K^{*}.
$$

Empirically, the normal matrix $G$ can be compressed into an HSS matrix.

In the implementation, a BF approximation $\widetilde{K}$ is first constructed, and fast matrix-vector products with $\widetilde{K}^{*} \widetilde{K}$ are used as black-box operations to construct an HSS approximation $\widetilde{G}$ of $G$.

The HSS matrix $\widetilde{G}$ is then factorized using the ULV factorization, which enables fast application of the inverse operator $\widetilde{F} = \widetilde{G}^{-1}$ through structured solves, without explicitly forming $\widetilde{G}^{-1}$. Therefore,

$$
K^{-1} = G^{-1} K^{*} \approx \widetilde{F} \widetilde{K}^{*}.
$$

The algorithm consists of three steps:

1. Construct a BF approximation $\widetilde{K}$ of $K$.
2. Construct an HSS approximation $\widetilde{G}$ of $G$ using the black-box construction algorithm, where the BF approximation $\widetilde{K}^{*} \widetilde{K}$ is used to perform fast matrix-vector products.
3. Factorize $\widetilde{G}$ using the ULV factorization to enable fast application of $\widetilde{F} = \widetilde{G}^{-1}$ through structured solves.

Related paper:

- Yingzhou Li and Jingyu Liu, *Approximate inversion of discrete Fourier integral operators via hierarchically semiseparable matrices*, in preparation.

## Getting Started

### Requirements

- MATLAB;
- [FastBF](https://github.com/YingzhouLi/FastBF.m);

### Clone the Repository

```bash
git clone https://github.com/JingyuLiuMath/FIO.git
cd FIO
```

### Install FastBF

```bash
cd extern
git clone https://github.com/YingzhouLi/FastBF.m.git
cd ..
```

### Set the MATLAB Path

Start MATLAB from the repository root and run:

```matlab
fio_startup();
```

### Run the HSS-Based Inversion Tests

From the repository root, run the one-dimensional test with

```matlab
run("test/test_FIO1D_inv_HSS.m");
```

and the two-dimensional test with

```matlab
run("test/test_FIO2D_inv_HSS.m");
```

The tests have been successfully run on Windows and Linux with:

- MATLAB: `R2023b`;
- FIO: `1.0.0`.

## Interface of the HSS-Based Inversion

The one-dimensional workflow uses

```matlab
result_bf = run_FIO1D_inv_CG(...);
result = run_FIO1D_inv_HSS(result_bf, ...);
```

The two-dimensional workflow uses

```matlab
result_bf = run_FIO2D_inv_CG(...);
result = run_FIO2D_inv_HSS(result_bf, ...);
```

The `run_FIO1D_inv_CG` and `run_FIO2D_inv_CG` functions construct the BF, check its accuracy, and run the unpreconditioned CG baseline. The `run_FIO1D_inv_HSS` and `run_FIO2D_inv_HSS` functions construct and factorize the HSS approximation, perform a direct approximate solve, and run the HSS-preconditioned PCG.

Currently, the HSS-based inversion supports both the SMP and BHP representations of BF in one and two dimensions.

## Reproducibility

### Inverse Based on HSS Matrices

The HSS experiments are located in [`experiments/fio_inv_hss`](./experiments/fio_inv_hss) and cover:

- one-dimensional constant-amplitude FIOs;
- one-dimensional variable-amplitude FIOs;
- two-dimensional constant-amplitude FIOs;
- two-dimensional variable-amplitude FIOs.

Each experiment follows two stages. For example, the one-dimensional constant-amplitude experiment is run with

```matlab
cd experiments/fio_inv_hss/1d_const_amplitude
exp_1d_const_amplitude_cg
exp_1d_const_amplitude_indep
```

The parameter grids and data paths are defined in the corresponding `exp_*_settings.m` file. On Unix systems, the default settings use paths under `/scratch/jyliu/FIO`. Users running on another system should update
`data_path` before starting an experiment.

## License

This project is released under the [MIT License](./LICENSE).

## Developers

- [Jingyu Liu](https://jingyuliumath.github.io/)
- [Yingzhou Li](https://yingzhouli.com/)
