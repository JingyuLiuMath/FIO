# FIO

## Overview

MATLAB research code for Fourier integral operators (FIOs).

## Getting Started

### Clone the Repository

``` bash
git clone git@github.com:JingyuLiuMath/FIO.git
cd FIO
```

### Install Related BF Packagres

See [readme](./extern/readme.md).

### Set MATLAB Path

Start MATLAB from the repository root and run:

``` matlab
fio_startup();
```

### Running the Tests

From the repository root, run:

``` matlab
run("test/test_FIO1D_inv_HSS.m");
```

The test has been successfully run on both Windows and Linux with the following software versions:

- MATLAB: `R2023b`;
- FIO: `1.0.0`.

## Reproducibility

The `experiments` contains examples discussed in the experimental sections of related papers.

- `fio_inv_hss`
  - Approximate inversion of discrete FIOs via HSS matrices.
  - Related paper: TBA

## Developers

- [Jingyu Liu](https://jingyuliumath.github.io/)
- [Yingzhou Li](https://yingzhouli.com/)
