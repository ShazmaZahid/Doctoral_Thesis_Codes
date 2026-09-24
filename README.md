
# $\star$-Methods for Computing the Action of the Matrix Exponential

![MATLAB](https://img.shields.io/badge/Language-MATLAB-orange.svg)

## 📌 Overview

This repository contains the MATLAB source code for the numerical experiments
of my doctoral thesis on **star-methods for computing the action of the matrix
exponential**.


## 📂 Project Structure

Each publication has its own directory containing the star-methods, the driver
that reproduces its results, and the helper functions it needs.

```
Doctoral_Thesis_Codes/
├── Paper1_PAMM2024/     base star-methods for exp(A)v          (conference proceedings)
└── Paper2_LAA2026/      Schur-based star-methods for exp(tA)v  (journal)
```

## ⚙️ Prerequisites

* **MATLAB** (tested on R2023b), with the **Control System Toolbox** (for `dlyap`).
* **[Chebfun](https://www.chebfun.org/)** — used by `genCoeffMatrix.m`.
* The Schrodinger examples (`test_id` 8–10 in Paper 2) additionally require the **PDE Toolbox**.

### Comparison routines (third party)

The star-methods are compared against two established, freely available
routines. They are **not redistributed here**; download them from their
authors and add them to the MATLAB path:

| Routine | Source | Reference |
| :--- | :--- | :--- |
| **`expv`** | [Expokit](https://www.maths.uq.edu.au/expokit/) (R. B. Sidje) | R. B. Sidje, *Expokit: A software package for computing matrix exponentials*, ACM Trans. Math. Software 24 (1998), 130–156. |
| **`expmv`, `expmv_tspan`** | [github.com/higham/expmv](https://github.com/higham/expmv) (Al-Mohy & Higham) | A. H. Al-Mohy & N. J. Higham, *Computing the action of the matrix exponential, with an application to exponential integrators*, SIAM J. Sci. Comput. 33 (2011), 488–511. |

`expmv_tspan` requires the full `higham/expmv` package (`expmv.m`,
`select_taylor_degree.m`, `normAm.m`, and the `theta_taylor*.mat` data files).

The drivers write their timings and relative errors to a dated text file,
`data_comp_time<date>.txt`.

## 📜 File Descriptions — Paper 1 (conference proceedings)

Pozza, S. & Zahid, S., *A new Legendre polynomial approach for computing the
matrix exponential action on a vector*, **Proceedings in Applied Mathematics
and Mechanics** 24(4) (2024), e202400049.
[DOI: 10.1002/pamm.202400049](https://doi.org/10.1002/pamm.202400049)

Reproduces **Tables 1 and 2** of the paper (Section 4.2 of the thesis).

| File Name | Description |
| :--- | :--- |
| **`first_star_method.m`** | First star-method. 
| **`second_star_method.m`** | Second star-method. 
| **`table1_first_star.m`** | Reproduces the whole of Table 1 in one run. |
| **`table2_second_star.m`** | Reproduces the whole of Table 2 in one run. |
| **`Arnoldimgs.m`** | Arnoldi iteration with modified Gram-Schmidt orthogonalization. |
| **`genCoeffMatrix.m`** | Coefficient matrix of $f(t)\Theta(t-s)$ in a basis of orthonormal Legendre polynomials. Requires chebfun. |
| **`California.mat`** | Adjacency matrix of Kleinberg's "California" web-search network (Pajek collection). |

## 📜 File Descriptions — Paper 2 (journal)

Pozza, S. & Zahid, S., *Computing the action of a matrix exponential on an
interval via the star-product approach*, **Linear Algebra and its
Applications** (2026), in press.
[DOI: 10.1016/j.laa.2026.04.014](https://doi.org/10.1016/j.laa.2026.04.014)

Reproduces **Tables 4.3–4.5** of the thesis (Section 4.3): 

| File Name | Description |
| :--- | :--- |
| **`experiments_paper2.m`** | Main driver. Selects one of ten benchmark problems (`test_id`) and reports, for each method, the relative error and the average computational time. |
| **`run_star_method.m`** | Base star-method|
| **`run_star_dlyap.m`** | Schur-based star-method (direct variant) |
| **`run_star_krylov.m`** | Schur-based star-method with Arnoldi |
| **`setup_matrix_and_vector.m`** | Builds the test matrix and vector for each of the ten benchmark problems. |
| **`Arnoldimgs.m`** | Arnoldi iteration with modified Gram-Schmidt orthogonalization. |
| **`genCoeffMatrix.m`** | Coefficient matrix of $f(t)\Theta(t-s)$ in a basis of orthonormal Legendre polynomials. Requires chebfun. |

## 🚀 Usage

1.  Install MATLAB and add the **Chebfun** library to the path.
2.  For Paper 2, download **`expv`** and the **`higham/expmv`** package (see the
    table above) and add them to the path.
3.  Clone this repository:
    ```bash
    git clone https://github.com/ShazmaZahid/Doctoral_Thesis_Codes.git
    ```
4.  Run a script from the relevant directory:
    ```matlab
    cd Paper1_PAMM2024
    table2_second_star        % reproduces Table 2

    cd ../Paper2_LAA2026
    experiments_paper2        % set test_id (1-10) at the top
    ```

The random-matrix examples use a fixed seed (`rng(0)`) where applicable, so
their results are reproducible; the deterministic examples reproduce the
published values.

## ✍️ Authors

* **Shazma Zahid** — *Doctoral Thesis*
* **Stefano Pozza, Dr., Ph.D.** — *Supervisor*
