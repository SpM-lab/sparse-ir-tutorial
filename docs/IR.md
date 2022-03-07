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

# Intermediate representation (IR)

$$
\newcommand{\iv}{{\mathrm{i}\nu}}
\newcommand{\wmax}{{\omega_\mathrm{max}}}
\newcommand{\dd}{{\mathrm{d}}}
$$

+++

## One-particle Green's function

We introduce a Green’s function with imaginary arguments in time and frequency.
This has no physical meaning but is merely a mathematical trick and makes calculations easier (to give another example of this: in Minsowski spacetime we take advantage of a similiar substitution)

The so-called imaginary-frequency (Matsubara) Green's functions are defined as followed: 

$$
G_{ij}(\tau-\tau')=-\langle T_\tau [c_i(\tau){c}^\dagger_j(\tau')]\rangle,
$$

where $i$ and $j$ denote spin/orbital/band and $T_\tau$ is the time-ordering operator.
Here $\tau$ represents a imaginary time unit $\mathrm{i}t$,
while $c_i$ and $c_j$ are the annihilation and creation operators.

The Fourier Transformation of $G_{ij}(\tau)$ (with $\tau \in [0,\beta]$) reads

$$
G_{ij}(\iv_n) = \int_0^{\beta} \dd \tau e^{\iv_n\tau} G_{ij}(\tau),
$$

where $\nu_n = (2n+1)\pi/\beta$ (fermion) and $\nu_n = 2n\pi/\beta$ (boson) with $n$ being an integer.
The inverse temperature is denoted by $\beta$ (We take $\hbar=1$).
The inverse transmation is given by

$$
G_{ij}(\tau) = \frac{1}{\beta}\sum_{n=-\infty}^\infty e^{-\iv_n\tau}G_{ij}(\iv_n).
$$ 

Continuing $G_{ij}(\iv_n)$ to a holomorphic function in the upper half of the complex plane,
the imaginary-frequency (Matsubara) Green's function can be related to the "orginary" retarded Green's function as

$$
G_{ij}^\mathrm{R}(\omega)=G_{ij}(z \rightarrow \omega+\mathrm{i}0^{+})
$$

In the following, we omit the symbols $i$, $j$, $n$ unless there is confusion.

+++

## Lehmann representation

In the imaginary-frequency domain, the Lehmann representation reads

$$
\begin{align}
    G(\iv) &= \int_{-\infty}^\infty \dd\omega \underbrace{\frac{1}{\iv - \omega}}_{\equiv K(\iv, \omega)} A(\omega),
\end{align}
$$

where $A(\omega)$ is a spectral function and we take $\hbar=1$.
$K(\iv, \omega)$ is the so-called analytic continuation kernel.
The Lehmann representation can be transformed to the imaginary-time domain as

$$
\begin{align}
    G(\tau) &= - \int_{-\infty}^\infty \dd\omega K(\tau, \omega) A(\omega),
\end{align}
$$ (lehmann-tau)

where $0 < \tau < \beta$ and 

$$
\begin{align}
    K(\tau, \omega) &\equiv - \frac{1}{\beta} \sum_{\iv} e^{-\iv \tau} K(\iv, \omega) =
    \begin{cases}
        \frac{e^{-\tau\omega}}{1+e^{-\beta\omega}} & (\mathrm{fermion}),\\
        \frac{e^{-\tau\omega}}{1-e^{-\beta\omega}} & (\mathrm{boson})
    \end{cases}.
\end{align}
$$

The minus sign originates from our the convention $K(\tau, \omega) > 0$.
To avoid the divergence of the bosonic kernel at $\omega=0$, we reformulate Eq. {eq}`lehmann-tau` as

$$
\begin{equation}
    G(\tau)= - \int_{-\infty}^\infty\dd{\omega} K^\mathrm{L}(\tau,\omega) \rho(\omega),
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
{cite:t}`Shinaoka:2017ix`

The singular value expnasion of the kernel {eq}`KL` reads

$$
\begin{align}
    K^\mathrm{L}(\tau, \omega) &= \sum_{l=0}^\infty U_l(\tau) S_l V_l(\omega),
\end{align}
$$

for $\omega \in [-\wmax, \wmax]$ with $\wmax$ (>0) being a cut-off frequency.

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
