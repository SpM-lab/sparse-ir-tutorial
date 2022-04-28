---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Sparse sampling

In this page, we describe how to infer IR expansion coefficients using the sparse-sampling techiniques.

```{code-cell} ipython3
import sparse_ir
import numpy as np
%matplotlib inline
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 15
```

## Setup
We consider a semicircular spectral modeli (full bandwidth of 2):

$$
\rho(\omega) = \frac{2}{\pi}\sqrt{1-\omega^2}.
$$

First, we compute the numerically exact expansion coefficients $g_l$.
Below, we plot the data for even $l$.

```{code-cell} ipython3
def rho(omega):
    if np.abs(omega) < 1:
        return (2/np.pi) * np.sqrt(1-omega**2)
    else:
        return 0.0

beta = 10000
wmax = 1
basis = sparse_ir.FiniteTempBasis("F", beta, wmax, eps=1e-15)

rhol = basis.v.overlap(rho)
gl = - basis.s * rhol

ls = np.arange(basis.size)
plt.semilogy(ls[::2], np.abs(gl[::2]), marker="x", label=r"$|g_l|$")
plt.xlabel(r"$l$")
plt.ylabel(r"$|g_l|$")
plt.ylim([1e-15, 10])
plt.legend()
plt.show()
```

## From sampling times

We first create a `TauSampling` object for the default sampling times.

```{code-cell} ipython3
smpl_tau = sparse_ir.TauSampling(basis)
print("sampling times: ", smpl_tau.sampling_points)
print("Condition number: ", smpl_tau.cond)
```

The condition number is around 50, indicating that 1--2 significant digits may be lost in a fit from the sampling times. Let us fit from the sampling times!

```{code-cell} ipython3
# Evaluate G(τ) on the sampling times
gtau_smpl = smpl_tau.evaluate(gl)
plt.plot(smpl_tau.sampling_points, gtau_smpl, marker="x", ls="")
plt.xlabel(r"$\tau$")
plt.ylabel(r"$G(\tau)$")
plt.show()

# Fit G(τ) on the sampling times
gl_reconst_from_tau = smpl_tau.fit(gtau_smpl)
```

## From sampling frequencies

We create a `MatsubaraSampling` object for the default sampling frequencies.

```{code-cell} ipython3
smpl_matsu = sparse_ir.MatsubaraSampling(basis)
print("sampling frequencies: ", smpl_matsu.sampling_points)
print("Condition number: ", smpl_matsu.cond)
```

The condition number is slightly larger than that for the sampling times.

```{code-cell} ipython3
# Evaluate G(iv) on the sampling frequencies
giv_smpl = smpl_matsu.evaluate(gl)
plt.plot((np.pi/beta)*smpl_matsu.wn, giv_smpl.imag, marker="x", ls="")
plt.xlabel(r"$\nu$")
plt.ylabel(r"Im $G(\mathrm{i}\nu)$")
plt.show()

# Fit G(τ) on the sampling times
gl_reconst_from_matsu = smpl_matsu.fit(giv_smpl)
```

## Comparison with exact results
We now compare the reconstructed expansion coefficients with the exact one. For clarity, we plot only the data for even $l$.

```{code-cell} ipython3
plt.semilogy(ls[::2], np.abs(gl[::2]), marker="", ls="-", label="Exact")
plt.semilogy(ls[::2], np.abs(gl_reconst_from_tau[::2]), marker="x", label="from sampling times")
plt.semilogy(ls[::2], np.abs(gl_reconst_from_matsu[::2]), marker="+", label="from sampling frequencies")
plt.xlabel(r"$l$")
plt.ylabel(r"$|g_l|$")
plt.ylim([1e-17, 10])
plt.legend(frameon=False)
plt.show()
```

We saw a perfect match! Let us plot the differences from the exact one. 

```{code-cell} ipython3
plt.semilogy(ls[::2], np.abs((gl_reconst_from_tau-gl)[::2]), marker="x", label="from sampling times")
plt.semilogy(ls[::2], np.abs((gl_reconst_from_matsu-gl)[::2]), marker="+", label="from sampling frequencies")
plt.xlabel(r"$l$")
plt.ylabel(r"Error in $g_l$")
plt.ylim([1e-18, 10])
plt.legend()
plt.show()
```
