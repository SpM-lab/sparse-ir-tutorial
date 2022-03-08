---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.7
kernelspec:
  display_name: Python 3.9.10 64-bit
  language: python
  name: python3910jvsc74a57bd0b0fa6594d8f4cbf19f97940f81e996739fb7646882a419484c72d19e05852a7e
---

+++ {"tags": []}

# Transformation from/to IR

$$
\newcommand{\iv}{{\mathrm{i}\nu}}
\newcommand{\wmax}{{\omega_\mathrm{max}}}
\newcommand{\dd}{{\mathrm{d}}}
$$

In this section, we explain how to transform numerical data to IR.

+++

## From smooth spectral function

It is straightforwad to compute expansion coefficients in the IR from a given spectral function $\rho(\omega)$.
Note that $\rho(\omega)$ is the modified spectral function for bosons.
The expansion coefficients can be evaluated by computing the integral

$$
\begin{align}
    \rho_l &= \int_{-\wmax}^\wmax \dd \omega V_l(\omega) \rho(\omega).
\end{align}
$$

One might consider to use the Gauss-Legendre quadrature.
As seen in previous sections, the distribution of $V_l(\omega)$ is much denser than Legendre polynomial $P_l(x(\tau))$ around $\tau=0, \beta$.
Thus, evaluating the integral precisely requires the use of composite Gauss–Legendre quadrature,
where the whole inteval $[-\wmax, \wmax]$ is divided to subintervals and the normal Gauss-Legendre quadrature is 
applied to each interval.
The roots of $V_l(\omega)$ for the highest $l$ used in the expansion
is a reasonable choice of the division points.
If $\rho(\omega)$ is smooth enough within each subinterval,
the result converges exponentially with increasing the degree of the Gauss-Legendre quadrature.

Below, we demonstrate how to compute $\rho_l$ for a spectral function consisting of of three Gausssian peaks using the composite Gauss-Legendre quadrature.
Then, $\rho_l$ can be transformed to $g_l$ by multiplying it with $- S_l$.

```{code-cell} ipython3
import sparse_ir
import numpy as np
import matplotlib.pyplot as plt

# Three Gaussian peaks (normalized to 1)
gaussian = lambda x, mu, sigma:\
    np.exp(-((x-mu)/sigma)**2)/(np.sqrt(np.pi)*sigma)

rho = lambda omega: 0.2*gaussian(omega, 0.0, 0.15) + \
    0.4*gaussian(omega, 1.0, 0.8) + 0.4*gaussian(omega, -1.0, 0.8)

omegas = np.linspace(-5, 5, 1000)
plt.xlabel(r"$\omega$")
plt.ylabel(r"$\rho(\omega)$")
plt.plot(omegas, rho(omegas))
plt.show()
```

```{code-cell} ipython3
beta = 10
wmax = 10
basis = sparse_ir.FiniteTempBasis("F", beta, wmax, eps=1e-10)

rhol = basis.v.overlap(rho)
gl = - basis.s * rhol

plt.semilogy(np.abs(rhol), marker="o", label=r"$|\rho_l|$")
plt.semilogy(np.abs(gl), marker="x", label=r"$|g_l|$")
plt.xlabel(r"$l$")
plt.ylim([1e-5, 1])
plt.legend()
plt.show()
```

$\rho_l$ is evaluated on arbitrary real frequencies as follows.

```{code-cell} ipython3
rho_omgea_reconst = basis.v(omegas).T @ rhol

plt.xlabel(r"$\omega$")
plt.ylabel(r"$\rho(\omega)$")
plt.plot(omegas, rho(omegas))
plt.show()
```

## From IR to imaginary time

We are now ready evaluate $g_l$ on arbitrary $\tau$ points.
A naive way is as follows.

```{code-cell} ipython3
taus = np.linspace(0, beta, 1000)
gtau1 = basis.u(taus).T @ gl
plt.plot(taus, gtau1)
plt.xlabel(r"$\tau$")
plt.ylabel(r"$G(\tau)$")
plt.show()
```

Alternatively, we can use the ``sampling`` module as follows.

```{code-cell} ipython3
smpl = sparse_ir.TauSampling(basis, taus)
gtau2 = smpl.evaluate(gl)
plt.plot(taus, gtau1)
plt.xlabel(r"$\tau$")
plt.ylabel(r"$G(\tau)$")
plt.show()
```

## From full imaginary-time data

A numerically stable way to compute the expansion coefficients of $G(\tau)$ in IR 
is evaluating the integral

$$
G_l = \int_0^\beta \dd \tau G(\tau) U_l(\tau).
$$

You can use `overlap` function as well.

```{code-cell} ipython3
def eval_gtau(taus):
    uval = basis.u(taus) #(nl, ntau)
    return uval.T @ gl

gl_reconst = basis.u.overlap(eval_gtau)

plt.semilogy(np.abs(gl_reconst), label="reconstructed", marker="o")
plt.semilogy(np.abs(gl), label="exact", marker="x")
plt.semilogy(np.abs(gl_reconst - gl), label="error", marker="")
plt.xlabel(r"$l$")
plt.xlabel(r"$|g_l|$")
plt.ylim([1e-20, 1])
plt.legend()
plt.show()
```

```{code-cell} ipython3

```

$$
G_l = \int_0^\beta \dd \tau G(\tau) U_l(\tau).
$$

We demonstrate this for $G(\tau)$ of the single Hubbard atom:

$$
G(\tau) = - \frac{1}{2}\left( \frac{e^{-\tau U/2}}{1 + e^{-\beta U/2}} + \frac{e^{\tau U/2}}{1 + e^{\beta U/2}} \right)
$$

with $\beta=300$ and $U=1$.

```{code-cell} ipython3
gtau_single_pole = lambda tau, epsilon: -np.exp(-tau*epsilon)/(1+np.exp(-beta*epsilon))
gtau = lambda taus: 0.5*(gtau_single_pole(taus, 0.5*U) + gtau_single_pole(taus, -0.5*U))
```

```{code-cell} ipython3

```
