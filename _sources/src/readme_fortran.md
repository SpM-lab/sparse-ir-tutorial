# sparse-ir-fortran - A Fortran interface for sparse-ir

This library provides a Fortran interface for `sparse-ir`.
At runtime, the installation of `sparse-ir` is not required.

First, please install `sparse-ir` with `xprec`.

## Creating a data file
### 1. Generate a data file:

```bash
> python3 dump.py 1e+4 1e-10 ir_nlambda4_ndigit10.dat
```

This generate the data file `ir_nlambda4_ndigit10.dat` containing sparse sampling points and transformation matrices for $\Lambda=10^4$ and $\epsilon = 10^{-10}$ ($\epsilon$ is a cut-off value for singular values).

### 2. Build object files and link them to your program:

```bash
> gfortran  -c -o sparse_ir.o sparse_ir.f90
> gfortran  -c -o sparse_ir_io.o sparse_ir_io.f90
```

You do not need `sparse_ir_preset.f90` if your program reads a data file at runtime.

## Embedding data in a Fortran source file
Data can be embedded in a source file.
This allows you to avoid loading data files at runtime.

### 1. Generate a fortran source file:

The following command generates a source file containg data for a matrix of
$\Lambda=10^{\mathrm{nlambda}}$ (nlambda = 1, 2, 3, 4) and $\epsilon=10^{-\mathrm{ndigit}}$ (ndigit = 10).

```bash
> python3 mk_preset.py --nlambda 1 2 3 4 --ndigit 10 > sparse_ir_preset.f90
```

### 2. Build object files and link them to your program:

The generated file `sparse_ir_preset.f90` is consisitent with Fortran95 and can be compiled as follows:

```bash
> gfortran  -c -o sparse_ir.o sparse_ir.f90
> gfortran  -c -o sparse_ir_preset.o sparse_ir_preset.f90
```

You do not need `sparse_ir_io.f90` if your program uses embedded data.

## How to use sparse-ir modules in your Fortran program
If your program reads data from the file, you should declare the use of the sparse-ir modules to initiate IR basis objects as follows:
```fortran
program main
  use sparse_ir
  use sparse_ir_io
```
Or, if you want to use the embed data of IR basis objects, you should do this instead:
```fortran
program main
  use sparse_ir
  use sparse_ir_preset
```
In either case, the derived type of "`IR`" should be declared.
```fortran
  implicit none
  type(IR) :: ir_obj ! You can declare with any name you wish.
```

Hereafter it is assumed that you will set the parameters as $\Lambda = 10^4$, $\beta = 10^3$, and $\epsilon = 10^{-10}$.
You can store the IR basis objects into the derived type of "`IR`" as follows:

Using `sparse_ir_io`:
```fortran
  double precision :: beta ! inverse temperature
  ...
  beta = 1.0d3
  open(99, file="ir_nlambda4_ndigit10.dat", status='old') ! Any unit number is OK.
  ir_obj = read_ir(99, beta)
```
Using `sparse_ir_preset`:
```fortran
  double precision :: beta ! inverse temperature
  ...
  beta = 1.0d3
  ir_obj = mk_ir_preset(4, 10, beta)
```

Here you are ready to use the IR basis objects and call the IR basis subroutines for a given value of `beta` in your program. Note that a derived type of `IR` (`ir_obj` here) should be updated by using the functions `read_ir` or `mk_ir_preset` each time you change the value of `beta`.

## Available objects

### `DOUBLE PRECISION:: IR%beta`
It returns is the input value $\beta$ of the functions `read_ir` or `mk_ir_preset`. <!--The type is `DOUBLE PRECISION`.-->
### `DOUBLE PRECISION:: IR%s (IR%size)`
It returns the singular values obtained from SVD of the kernel matrix. <!--The type is 1-dimensional `DOUBLE PRECISION` array and its size is `IR%size`.-->
### `INTEGER:: IR%size`
It returns the size of `IR%s`. <!--The type is integer.-->
### `DOUBLE PRECISION:: IR%tau (IR%ntau)`
It returns the values of $\tau$ of sparse sampling points. <!--The type is 1-dimensional array and its size is `IR%ntau`.-->
### `INTEGER:: IR%ntau`
It returns the size of `IR%tau`. <!--The type is integer.-->
### `INTEGER:: IR%freq_f (IR%nfreq_f)`
It returns the odd integers corresponding to sampling Matsubara frequencies for fermionic function. <!--The type is 1-dimensional array and its size is `IR%nfreq_f`.-->
### `INTEGER:: IR%nfreq_f`
It is the number of sampling Matsubara frequencies for fermionic function.
### `INTEGER:: IR%freq_b (IR%nfreq_b)`
It returns the even integers corresponding to sampling Matsubara frequencies for bosonic function. <!--The type is 1-dimensional array and its size is `IR%nfreq_b`.-->
### `INTEGER:: IR%nfreq_b`
It is the number of sampling Matsubara frequencies for bosonic function.

### `TYPE(DecomposedMatrix):: IR%uhat_f`
It refers the derived type of `DecomposedMatrix` which contains `IR%uhat_f%a`, `IR%uhat_f%inv_s`, `IR%uhat_f%ut`, and `IR%uhat_f%v`. When TYPE "`IR`" is defined for a given `beta`, SVD of $\{\hat{U}_l(\mathrm{i}\nu_n)\}$ for the `beta` is performed to define `IR%uhat_f%inv_s`, `IR%uhat_f%ut`, and `IR%uhat_f%v`, which are used in subroutines `fit_matsubara_f` and `evaluate_matsubara_f`.

The basis functions on fermionic sampling Matsubara frequencies is SVDecomposed in advance  as follows:
$$
\begin{align*}
\hat{U}_l(\mathrm{i}\nu_n) &= A_{nl} \\
&=\sum_{r,r'}U_{nr} \Sigma_{rr'} (V^\mathrm{T})_{r'l} \\
(A &= U \Sigma V^\mathrm{T}),
\end{align*}
$$
with
$$
\Sigma_{rr'} = \sigma_r\delta_{r,r'}~  (\sigma_r > \epsilon).
$$
If a following fitting problem is given for a certain fermionic function $G(\mathrm{i}\nu_n)$ defined on Matsubara frequencies,
$$
\sum_{l}\hat{U}_l(\mathrm{i}\nu_n)G_l\approx G(\mathrm{i}\nu_n),
$$
you can solve the problem using `fit_matsubara_f` as follows:
$$
G_l\approx \sum_{r, n}V_{lr}\Sigma^+_{rr}(U^\mathrm{T})_{rn}G(\mathrm{i}\nu_n).
$$
`IR%uhat_f%a`, `IR%uhat_f%ut`, and `IR%uhat_f%v` are 2-dimensional arrays the matrice corresponding to $A$, $U^\mathrm{T}$, and $V$, respectively. `IR%uhat_f%inv_s` is the 1-dimensional array storing the components of the diagonal matrix $\Sigma^+$, namely $1/\sigma_r$. `IR%uhat_f%ns` is the size of  `IR%uhat_f%inv_s`. `IR%uhat_f%m` and `IR%uhat_f%n` equal to `ir_obj%nfreq_f` and `ir_obj%size`, respectively.

<!---
### `IR%uhat_f%a` (`IR%uhat_b%a`)
It is the IR basis set on Matsubara frequencies, $\{\hat{U}_l(\mathrm{i}\nu_n)\}$. `IR%uhat_f%a` is for fermionic functions and `IR%uhat_b%a` is for bosonic functions. The type is 2-dimensional array and its shape is `(IR%uhat_f%m, IR%uhat_f%n)` which equals to `(IR%nfreq_f, IR%size)`.

### `IR%uhat_f%ns` (`IR%uhat_b%ns`)
It returns the number of size of  `IR%uhat_f%inv_s`.
-->

## Subroutines

### `SUBROUTINE fit_matsubara_f`
The subroutine fits a set of expansion coefficients <!--of the IR basis -->$G_l$ to a given fermionic function $G(\mathrm{i}\nu_n)$ on sampling Matsubara frequencies by using SVD.
$$
\begin{align*}
G_l = {\mathop{\rm argmin}\limits}_{G_l}\left|G(\mathrm{i}\nu_n) - \sum_{l}\hat{U}_l(\mathrm{i}\nu_n)G_l \right|^2
\end{align*}
$$
The Usage is
```fortran
call fit_matsubara_f(obj, g_in, g_out)
```
where `g_in` and `g_out` correspond to $G(\mathrm{i}\nu_n)$ and $G_l$, respectively. `obj` is the derived type of "`IR`", and `g_in` and `g_out` are 2-dimensional `COMPLEX(KIND(0D0))` arrays. The inputs are `obj` and `g_in` and the output is `g_out`. Before calling this subroutine, you should reshape the array of $G(\mathrm{i}\nu_n)$ to a 2-dimensional array whose last axis corresponds to $l$ and allocate `g_out` with appropriate shape.
That is, `g_in` and `g_out` should be allocated so as to have shapes of `(**, obj%nfreq_f)` and `(**, obj%size)`, respectively.

### `SUBROUTINE evaluate_matsubara_f`
This subroutine reconstructs the fermionic function $G(\mathrm{i}\nu_n)$ on sampling Matsubara frequencies from a given set of expansion coefficients $G_l$ as follows:
$$
\begin{align*}
G(\mathrm{i}\nu_n) = \sum_{l}\hat{U}_l(\mathrm{i}\nu_n)G_l
\end{align*}
$$
The Usage is
```fortran
call evaluate_matsubara_f(obj, g_in, g_out)
```
where `g_in` and `g_out` correspond to $G_l$ and $G(\mathrm{i}\nu_n)$, respectively. `obj` is the derived type of "`IR`", and `g_in` and `g_out` are 2-dimensional `COMPLEX(KIND(0D0))` arrays. The inputs are `obj` and `g_in` and the output is `g_out`. Before calling this subroutine, you should reshape the array of $G_l$ to a 2-dimensional array whose last axis corresponds to $\mathrm{i}\nu_n$ and allocate `g_out` with appropriate shape.
That is, `g_in` and `g_out` should be allocated so as to have shapes of `(**, obj%size)` and  `(**, obj%nfreq_f)`, respectively.