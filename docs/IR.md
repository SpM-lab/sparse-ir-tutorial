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

# Intermediate representation

$\newcommand{\wmax}{{\omega_\mathrm{max}}}$
$\newcommand{\dd}{{\mathrm{d}}}$

+++

## Lehmann representation

In the imaginary-time domain, the Lehmann represenation reads

$$
\begin{align}
    G(\tau) &= - \int_{-\infty}^\infty \dd\omega K(\tau, \omega) A(\omega),
\end{align}
$$

where $0 < \tau < \beta$, $A(\omega)$ is a spectral function, $K(\tau, \omega)$ is a *kernel*.
We take $\hbar = 1$.

The kernel is defined as

$$
\begin{align}
    K(\tau, \omega) &=
    \begin{cases}
        \frac{e^{-\tau\omega}}{1+e^{-\beta\omega}} & (\mathrm{fermion}),\\
        \frac{e^{-\tau\omega}}{1-e^{-\beta\omega}} & (\mathrm{boson})
    \end{cases}.
\end{align}
$$

To avoid the divergence of the bosonic kernel at $\omega=0$, we reformulate the Lehmann representation as

$$
\begin{equation}
    G(\tau)= - \int_{-\infty}^\infty\dd{\omega'} K^\mathrm{L}(\tau,\omega') \rho(\omega'),\label{eq:friedholm2}
\end{equation}
$$

where $K^\mathrm{L}(\tau, \omega)$ is the logistic kernel defined as

$$
K^\mathrm{L}(\tau, \omega) =  \frac{e^{-\tau\omega}}{1+e^{-\beta\omega}},
$$ (KL)

and $\rho(\omega)$ is the modified spectral function

$$
\begin{align}
    \rho(\omega) &\equiv 
    \begin{cases}
        A(\omega) & (\mathrm{fermion}),\\
        \frac{A(\omega)}{\tanh(\beta \omega/2)} & (\mathrm{boson}).
    \end{cases}
\end{align}
$$

+++

## Singular value expansion

{cite:p}`Shinaoka:2017ix`

{cite}`Shinaoka:2017ix`

The singular value expnasion of the kernel {eq}`KL` reads

$$
\begin{align}
    K^\mathrm{L}(\tau, \omega) &= \sum_{l=0}^\infty U_l(\tau) S_l V_l(\omega),
\end{align}
$$

for $\omega \in [-\wmax, \wmax]$.

+++

## Basis functions in imaginary frequency

Transformation to frequency

```{code-cell} ipython3
print(1)
```

```{code-cell} ipython3
print("AAAA")
```

```{code-cell} ipython3

```
