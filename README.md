# ⋆-Methods for Computing the Action of the Matrix Exponential

![MATLAB](https://img.shields.io/badge/Language-MATLAB-orange.svg)
![License](https://img.shields.io/badge/License-MIT-blue.svg)

## 📌 Overview

This repository contains the source code for the numerical experiments of my
doctoral thesis on **star-methods for computing the action of the matrix
exponential**.

The methods recast the computation of the matrix exponential action on a
vector as the solution of a Stein matrix equation, obtained through a
Legendre-polynomial discretization of the so-called star-product framework.
The code is organized by the publication or thesis chapter it reproduces.

## 📂 Project Structure

The implementation is written in **MATLAB**. Each publication has its own
directory containing the methods, the scripts that reproduce its tables, and
the helper functions it needs.

```
Doctoral_Thesis_Codes/
└── Paper1_PAMM2024/     base star-methods for exp(A)v
```

*(Codes for the second paper and for the preconditioner chapter will be added
in their own directories.)*

### Prerequisites

For the optimal execution of this code, the following are required:
* **MATLAB** (tested on R2023b), with the **Control System Toolbox** (for `dlyap`).
* **[Chebfun](https://www.chebfun.org/)** – An open-source software system for numerical computing with functions, used by `genCoeffMatrix.m`.

## 📜 File Descriptions — Paper 1 (PAMM 24, 2024)

Pozza, S. & Zahid, S., *A new Legendre polynomial approach for computing the
matrix exponential action on a vector*, **Proc. Appl. Math. Mech.** 24 (2024),
e202400049. [DOI: 10.1002/pamm.202400049](https://doi.org/10.1002/pamm.202400049)

The `Paper1_PAMM2024` directory contains the following scripts and functions:

| File Name | Description |
| :--- | :--- |
| **`first_star_method.m`** | First star-method (Section 4.1). Projects `A` onto a Krylov subspace by the Arnoldi process and solves the reduced Stein equation with `dlyap`. Prints the relative errors of the classical Krylov approximation and of the first star-method, as in Table 1. |
| **`second_star_method.m`** | Second star-method (Section 4.2). Solves the full star-linear system iteratively with GMRES, using a matrix-free operator so that the Kronecker matrix is never assembled. Prints the relative error and iteration count and plots the GMRES residual, as in Table 2 and Figure 1. |
| **`table1_first_star.m`** | Runs every case of Table 1 in a single run and prints the full table. |
| **`table2_second_star.m`** | Runs every case of Table 2 in a single run and prints the full table. |
| **`Arnoldimgs.m`** | Arnoldi iteration with modified Gram-Schmidt orthogonalization. Returns the orthonormal Krylov basis and the upper-Hessenberg matrix. |
| **`genCoeffMatrix.m`** | Computes the coefficient matrix of $f(t)\Theta(t-s)$ (where $\Theta$ is the Heaviside theta function) within a basis of orthonormal Legendre polynomials. |
| **`California.mat`** | Adjacency matrix of Kleinberg's "California" web-search network (Pajek collection), used in the California example. |

## 🚀 Usage

1.  Ensure you have MATLAB installed.
2.  Install the **Chebfun** library and add it to your MATLAB path.
3.  Clone this repository:
    ```bash
    git clone https://github.com/ShazmaZahid/Doctoral_Thesis_Codes.git
    ```
4.  Navigate to the `Paper1_PAMM2024` directory.
5.  Run a `table*.m` script to reproduce a full table, or run a single-method
    script after uncommenting one of its example blocks (Random, Poisson,
    California) at the top.

    ```matlab
    cd Paper1_PAMM2024
    table2_second_star        % prints Table 2
    ```

The random-matrix examples use a fixed seed (`rng(0)`) in the `table*.m`
scripts so their results are reproducible; the Poisson and California examples
are deterministic and reproduce the published values.

## ✍️ Authors

* **[Shazma Zahid]** – *Doctoral Thesis*
* **[Stefano Pozza, Dr., Ph.D.]** – *Supervisor*
