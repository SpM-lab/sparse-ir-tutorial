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

# GW approximation

+++

## Single Particle Green's function

+++

In order to give a short introduction into Green's functions in the context of many-body physics and in particular the GW approximation we will take a look at the single particle Green's function. It offers a relatively easy to understand solution for problems similar to the one we will deal with within the GW approximation and is thus very helpful in understanding the concept and idea behind such approximations. 

Starting off with a standard problem of quantum mechanics, the propagation of a single particle, the time-dependant Schrödinger's Equation gives the well-known formula for a particles wave function.

$$[i\partial_t-\mathcal{H}_0(r)-V(r)]\Psi(r,t)=0  $$  

It is now possible to define two (Green’s) functions, to which the following applies. 

$$\mathcal{H}_0(r)G_0(r,r’;t,t’)=\delta(r-r’)\delta(t-t’)$$<br>
$$[i\partial_t-\mathcal{H}_0(r)-V(r)]G_0(r,r';t,t')=\delta(r-r')\delta(t-t')$$

As we will see soon this will allow to obtain a solution for both a non-interacting particle ($G_0$) and a particle of some pertubation V, which is interacting with the system ($G$). Taking this further and making use of the inverse of the Green’s functions the solution for $\Psi$ can easily found to be 

$$\Psi(r,t)=\Psi^0(r,t)+\int dr' \int dt' G(r,r';t,t')V(r')\Psi^0(r',t')$$

for the full (interacting) Green's function and 

$$\Psi(r,t)=\Psi^0(r,t)+\int dr' \int dt' G_0(r,r';t,t')V(r')\Psi^0(r',t')$$

for a non-interacting particle. 

After iterating this solution for the wave function, it will lead to the so-called Dyson Equation with which we can write the full Green’s function $G$ as<br>

$$G=G_0 +G_0VG $$

[<sup id="fn1-back">1</sup>](#fn1).The Dyson Equation will be used later on to connect the Green's function with the self-energy $\Sigma$.

If we interpret the Green’s function as a propagator of the wave function $\Psi$, which is rather intuitiv since 

$$\Psi(r,t)=\int dr' G(r,r',t,t')\Psi(r',t')=\int dr' \langle r|e^{-i\mathcal{H}}(t-t')|r\rangle \langle r'|\Psi| r \rangle$$

We can define a function:

$$G(r,r',t-t’)= -i\theta (t-t')\langle r |e^ {-iH(t-t')} |r'\rangle$$

$G(r,r‘,t-t‘)$, often called retarded Green’s function[<sup id="fn2-back">2</sup>](#fn2), is nothing else than the probability amplitude for a particle to enter state r at a time t if it was injected into the system in the state r’ at time t’. This becomes even more evident if we take notice of the fact that in the Dyson equation we are basically calculating $G$ as the sum of all possible paths a particle can take in the system with $G_0$ as a non-interacting propagator and $V$ as a scattering event of some sort.[1][3]


There are of course higher order Green's functions to describe multi-particle-processes, these however will not be further explained as it would go beyond the scope.
 
[<sup id="fn1">1</sup>](#fn1-back) The solution is being obtained from $G=G_0 + G_0VG_0+ G_0VG_0VG_0+...=G_0 + G_0V(G_0+G_0VG)$.

[<sup id="fn2">2</sup>](#fn2-back)There is also a so-called advanced Green's function with $-i\theta(t'-t)$, which is used to describe the propagation of a quasiparticle or hole.

+++

## Matsubara Green's function

+++

We now introduce a Green’s function with imaginary arguments in time and frequency. This has no physical meaning but is merely a mathematical trick and makes calculations easier (to give another example of this: in Minsowski spacetime we take advantage of a similiar substitution)

The so called Matsubara Green's functions are defined as followed: 

$$G_{ij}(\tau-\tau')=-\langle T_\tau [c_i(\tau){c}_j(\tau')]\rangle$$

Here $\tau$ represents a imaginary time unit $it$, while $c_i$ and $c_j$ are the annihilation and creation operators.[<sup id="fn3-back">3</sup>](#fn3) With 

$$T_\tau [c_i(\tau)c_j(\tau')]= \theta(\tau-\tau')c_i(\tau)c_j(\tau')\pm\theta(\tau'-\tau)c_j(\tau')c_i(\tau)
\begin{cases}
+ \space for\space bosons\\
- \space for\space fermions\\
\end {cases}$$

the operators are arranged according to the course of events, which is why $T$ is often being refered to as the time ordering symbol in imaginary time. To be more precise: the left operator represents the later in time (Note that $T$ is differently defined for fermions and bosons). Similarly to what the single particle Green's function represents $G_{ij}(\tau-\tau')$ describes an probability amplitude to find a particle in the state $i$ at the time $\tau$ if it was in a state $j$ at time $\tau'$, at least for $\tau>\tau'$. If $\tau'>\tau$ applies we are again looking at the propagation of a quasi-particle (hole). 

Taking a look at the Fourier Transformation of $G_{ij}(\tau)$ (with $\tau \in [-\beta,\beta]$) it can be said that 

$$G_{ij}(i\omega_n)=\frac{1}{2}\int_{-\beta}^{\beta} d\tau  e^{i\omega_n\tau} G_{ij}(\tau)$$

and

$$G_{ij}(\tau) = \frac{1}{\beta}\sum_{n=-\infty}^\infty e^{-i\omega_n\tau}G_{ij}(i\omega_n).$$ 

As we can simplify this by separating the integral of $G_{ij}(i\omega_n)$ and making use of the fact that $n$ is odd or even depending on $(1\pm e^{-i\pi n})$ being $0$ or $2$, we can write 

$$G_{ij}(i\omega_n)=\int_0^{\beta} d\tau e^{i\omega_n\tau}G_{ij}(\tau)\space\space
\begin{cases}
\omega_n = \frac{2\pi n}{\beta}\space    for\space fermions \\
\omega_n =\frac{(2n+1)\pi}{\beta}\space  for\space bosons
\end{cases}$$

Since the functions $G_{ij}$ are referred to as Matsubara Green's functions, the discrete frequencies $\omega_n$ are called Matsubara frequencies. Because these frequencies contain $\beta$, which is essentially $\frac{1}{k_BT}$, information about the temperature is preserved.
            
The importance of the Matsubara Green's function in the context of many-body systems lies within its relation to the 'ordinary' retarded Green's function: If we imagine the complex plane, there exists a complex Green's function $G(z)$ whose imaginary part on the imaginary axis equals $G_{ij}(i\omega_n)$ (Matsubara Green's function) and whose real part equals a function $G_{ij}(\omega)$. As can be proven by the use of Lehmann's representation that function $G_{ij}(\omega)$ happens to be the retarded Green's function for the system. This makes it possible to compute $G_{ij}(i\omega_n)$ if the retarded Green's function is known and vice versa. That as well as analytic continuation allows us to state:

$$G_{ij}^R(\omega)=G_{ij}(i\omega_n\to\omega+i0^{+})$$

(for this $G_{ij}(z)$ has to be reduced to a rational function and be analytic, otherwise it will yield wrong results). An analog identity is true for the advanced Green's function.[1][2][3]

[<sup id="fn3">3</sup>](#fn3-back) $\langle A\rangle$ is the grand-canonical expectation value 
$$\langle A\rangle=\frac{1}{Z}T_r(e^{-\beta\mathcal{H}})$$

+++

Singular Value Expansion
------------------------

+++

Let us now take a closer look at $G_{ij}(i\omega_n)$. We can write our Matsubara Green's function as

$$\hat{G}^\alpha (i\omega)=\int_0^{\beta} d\tau e^{i\omega_n\tau}G(\tau)= \int_0^\beta d\tau e^{i\omega_n\tau} \langle T c_i(\tau)c_j(0)\rangle .$$

$\alpha$ denotes whether we are working with fermionic (F) or bosonic statistics (B). Another equivalent way of expressing $G_{ij}(i\omega_n)$ is through the sprectral Lehmann representation:

$$\hat{G}(i\omega) = \int_{-\omega_{max}}^{\omega_{max}} d\omega'\underbrace{ \frac{\omega'^{\delta_{\alpha,B}}}{i\omega-\omega'}}_{K^{\alpha}(i\omega,\omega')} \rho^\alpha (\omega')$$

Here $\rho^\alpha(\omega')$ represents the associated spectral density function in the real frequency domain $\omega_i \in [-\omega_{max},\omega_{max}]$ (we recall that these are connected to the imaginary frequencies $i\omega_n$ by analytic continuation). $\omega'^{\delta_{\alpha,B}}$ evaluates to $1$ or $\omega'$, depending on $\alpha=F$ or $\alpha=B$.

Because we are looking at an integral of the form $\int G(s,t)f(t)dt=u(s)$ (a so called Fredholm integral equation of the first kind) it is possible to perform a singular value expansion (SVE), so that $K^\alpha(i\omega,\omega')$ satisfies: 

<br>
$$K^\alpha(i\omega,\omega')=\sum_{l=0}^\infty \hat{U}_l^\alpha(i\omega)S_l^\alpha V(\omega')=\sum_{l=0}^{\infty}\hat{U}_l^\alpha(\tau)S_l^\alpha V_l^\alpha(\omega')$$


Here the representational form for both the imaginary-frequency and imaginary time domain can be seen.[1][4][7] In order to go further into detail we will construct that kernel $K$ (for both fermionic and bosonic statistics) using [sparse_ir](https://github.com/SpM-lab/sparse-ir) with the following conditions for our system:

```{code-cell}
import numpy as np
import matplotlib.pyplot as pl
import sparse_ir
import sys 
import math
```

```{code-cell}
T=0.1
wmax =1

beta = 1/T
```

```{code-cell}
#Fermionic Basis

basisf = sparse_ir.FiniteTempBasis('F', beta, wmax)

matsf=sparse_ir.MatsubaraSampling(basisf)
tausf=sparse_ir.TauSampling(basisf)
```

Intermediate Representation
------------------------

+++

We have now created a kernel $K^F(i\omega,\omega')$ with $U_l^F$ and $V_l^F$ as left and right singular functions and $S_l$ as the singular values (with $S_0>S_1>S_2>...>0$). The two  singular functions $U$ and $V$ make up the basis functions of the so-called Intermediate Representation (IR), which depends on $\beta$ and the cutoff $\omega_{max}$ as well as statistical properties (in this case fermionic properties $F$). Some of the information regarding real-frequency properties of the system is often lost during transition into the imaginary-time domain, so that the Matsubara Green's function does hold less information than the real-frequency Green's function. The reason for using IR lies within its compactness and ability to display that information in imaginary quantities.[5][7]

Additionally we calculated the positions of our sampling points for imaginary frequencies and imaginary time.

+++

The basis functions $V_l^F(\omega')$  (the right singular functions) form an orthonomal set on the real frequency axis.

```{code-cell}
wlins= np.linspace(-wmax, +wmax, 101)
pl.plot(wlins,basisf.v[0](wlins))
pl.plot(wlins,basisf.v[1](wlins))
pl.plot(wlins,basisf.v[18](wlins))
pl.legend(['l=0','l=1','l=18'])
pl.xlabel('$\omega$')
pl.ylabel('$V_l^F(\omega)$')
pl.ylim(-2,2)
pl.title('$V_l^F(\omega)$ over $\omega$ for different l')
```

On the other side we have $U_l^F(i\omega)$ which exists on the imaginary axis.

```{code-cell}
print('Matsubara frequencies :',matsf.wn)

pl.plot(matsf.wn,'+')
pl.title('Matsubara frequencies')
pl.xlabel('n')
pl.ylabel('$\omega_n$')
```

```{code-cell}
pl.plot(matsf.wn,basisf.uhat[0](matsf.wn).imag)
pl.plot(matsf.wn,basisf.uhat[1](matsf.wn).imag)
pl.plot(matsf.wn,basisf.uhat[18](matsf.wn).imag)
pl.legend(['l=0','l=1','l=18'])
pl.title('Imaginary part of $U_l^F(i\omega_n)$ over $i\omega_n$ for different l')
pl.xlabel('$\omega_n$')
pl.ylabel('$U_l^F(i\omega_n)$')
```

As we can express $K$ in both, the frequency and time domain, there is also a need for $U_l^F(\tau)$ (this will be further explained below)

```{code-cell}
print('Sampling points in tau:', tausf.tau)
pl.plot(tausf.tau,'x')
pl.title('Sampling points in tau')
```

```{code-cell}
pl.plot(tausf.tau,basisf.u[0](tausf.tau))
pl.plot(tausf.tau,basisf.u[1](tausf.tau))
pl.plot(tausf.tau,basisf.u[18](tausf.tau))
pl.legend(['l=0','l=1','l=18'])
pl.title('$U_l^F(tau)$')
```

In order to perform the GW approximation the same objects as above are also needed with regards to bosonic statistics.

```{code-cell}
#Bosonic Basis

basisb = sparse_ir.FiniteTempBasis('B', beta, wmax)
      
matsb=sparse_ir.MatsubaraSampling(basisb)
tausb=sparse_ir.TauSampling(basisb)
```

Interpreting $\hat G$ (with 'hat' denoting frequency representations) in terms of basis functions and expansion coefficients leads us to the following: 

Since $\hat{G}^\alpha(i\omega)$ has taken the form 

$$\hat G(i\omega)=\int K(i\omega,\omega')\rho(\omega')d\omega'$$

we can write

$$\hat G^\alpha(i\omega)=\sum_{l=0}^{N-1} \hat{U}_l^\alpha(i\omega_n^\alpha)\underbrace{S_l^\alpha \int V_l^\alpha(\omega')\rho(\omega')d\omega'}_{G_l^\alpha}=\sum_{l=0}^{N-1} \hat{U}^\alpha_l(i\omega_n^\alpha)G_l^\alpha$$.

Here we can see that $\hat{U}_l(i\omega_n)$ act as imaginary frequency basis functions with $G_l$ as the expansion coefficients. We can fourier transform $\hat{U}_l(i\omega)$  into imaginary time functions $U_l(\tau)$. These will simply be the basis functions of $G(\tau)$ in the imaginary time domain 

$$G^\alpha(\tau)=\sum_{l=0}^{N-1} U_l^\alpha(\bar\tau_k^\alpha) G_l^\alpha$$

This essentially means that once the IR expansion coefficients are known we are able to evaluate data at any frequency and time and can easily switch back and forth between the time and the frequency domain.[5]

In the following code we will assume $\rho(\omega)$ to be of the form:

$$\rho(\omega)=\frac{2}{\pi}*\sqrt{1-\omega^2}$$

with

$$\rho_l^F=\int_{-\omega_max}^{\omega_max}V_l^F(\omega')\rho(\omega')$$

```{code-cell}
def rho(x):
    return 2/np.pi*np.sqrt(1-x*x)

rho_l=basisf.v.overlap(rho)
#evaluate overlap \int dx p(x) f(x) using piecewise Gauss-Legendre quadrature
#p(x) are the polynomials
    
```

Once $\rho_l^F$ is known the expansion coefficients $G_l^F$ can be calculated using the singular values $S_l$.

$$G_l^F=S_l^F\rho_l^F $$

```{code-cell}
G_l_0=basisf.s*rho_l   
```

```{code-cell}
pl.semilogy(np.abs(basisf.s), ':')
pl.semilogy(np.abs(rho_l), '+-')
pl.semilogy(np.abs(G_l_0),'+-')
pl.legend(['basis S','rho_l','G_l'])
pl.title('Comparison of rho_l and $G_l$ over l')
pl.xticks(range(0,20))
pl.xlabel('l')
```

It can be seen that $G_l$ is getting smaller with growing $l$. This agrees with the fact that $S_{l=0}$ is the greatest singular value and $S_l$ gets smaller with growing $l$.

+++

We will now build $\hat G^F (i\bar{\omega}_n^F)$ using the expansion coefficients we just computed as well as the basis functions $\hat U_l^F(i\bar{\omega}_n^F)$ (we are still working in the fermionic basis).[5][7]


$$\hat G(i\bar{{\omega}}_n^F)=\sum_{l=0}^{N-1}G_l^F\hat U_l^F(i{\omega}_n^F)$$

```{code-cell}
#We compute G_iw two times as we will need G_iw_0 as a constant later on

G_iw_0=np.einsum('lv,l->v',basisf.uhat(matsf.wn),G_l_0)
G_iw_f=np.einsum('lv,l->v',basisf.uhat(matsf.wn),G_l_0)
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5))
ax[0].plot(matsf.wn,G_iw_f.imag)
ax[0].set_title('Imaginary part of $G(i\omega)$')
ax[1].plot(matsf.wn,G_iw_f.real)
ax[1].set_title('Real part of $G(i\omega)$')
```

As stated above transforming the Green's function from $\hat{G}$ into its time domain form  can be done rather easily. First we need to obtain $G_l$ by fitting it with e.g. the Least-Square method such that 

$$G_l^F=\underset{G_l}{\arg\min}   \sum_{i\bar{\omega}\in W^F}|\hat{G}(i\bar{{\omega}}_n^F)-\sum_{l=0}^{N-1}\hat{U}_l^F(i\bar{\omega}_n^F)G_l|^2$$

However, as we already computed $\hat G$ from the coefficients $G_l$, fitting is not really needed here but will be done for the sake of completeness.[5][7]

```{code-cell}
#Calculation of G_l_f using Least-Squares Fitting

G_l_f=matsf.fit(G_iw_f)
```

For the final step evaluation of $G_l^F$ on the (fermionic) $\tau$-sampling points is necessary.

$$G(\bar{\tau}_k^F)= \Sigma_{l=0}^{L-1}G_l^F U_l^F(\bar{\tau}_k^F)  $$

```{code-cell}
#G_tau,fermionic   

G_tau_f=np.einsum('lt,l->t',basisf.u(tausf.tau),G_l_f)
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[1].plot(tausf.tau,G_tau_f.imag,'+-')
ax[1].set_title('Imaginary part of $G(tau)$')
ax[0].plot(tausf.tau,G_tau_f.real,'+-')
ax[0].set_title('Real part of $G(tau)$')
```

The GF2 & GW method
-----------

+++

Above we have successfully derived $G^F(\bar{\tau}_k^F)$ from $\hat {G}^F(i\bar{\omega}_n^F)$ by first obtaining the Green's function's basis representation $G_l^F$ and secondly evaluating $G_l^F$ on the $\tau$-sampling points, which are given by the computations of sparse_ir. Here, the self-consistent second order Green's function theory (GF2) and the self-consistent GW methodare are based this very mechanism (and the so called Hedin equations).[5][7]

The following image should give an overview on their iterative schemes[5]:

<br>

![GF2_GW_pic.PNG](attachment:GF2_GW_pic.PNG)

+++

Evaluation and transformation in GW and GF2 is done in the same way (fitting and evaluating) for the quantities $G,P,W$ and $\Sigma$. The transition between these quantities is given by the Hedin equations. There are five Hedin equations in total. Because the vertex function (between $P$ and $W$) is approximated as $\Gamma=1$, we will only deal with four of them here.[9]

+++

As we already have computed $G_l^F$, the next step in order to perform a full iteration of GF2 & GW is the evaluation of the basis representation $G_l^F$ on the $\tau$-sampling points.  Contrary to before were we computed $G(\bar{\tau}_k^F)$, these sampling point will be the bosonic sampling points $\tau_k^B$. This allows us to switch to bosonic statistic and thus be able to calculate associated quantities like the Polarization $P$.[5]

$$G(\bar{\tau}_k^B)= \sum_{l=0}^{L-1}G_l^F  U_l^F(\bar{\tau}_k^B)$$

```{code-cell}
#G_tau, bosonic

G_tau_b=np.einsum('lt,l->t',basisb.u(tausb.tau),G_l_f)
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[1].plot(tausb.tau,G_tau_b.imag)
ax[1].set_title('Imaginary part of $G(tau)$')
ax[0].plot(tausb.tau,G_tau_b.real)
ax[0].set_title('Real part of $G(tau)$')
```

The so-called Polarization $P(\bar\tau_k^B)$ is given by the random phase approximation, a connection between two Green's functions:

$$P(\bar{\tau}_k^B)=G(\bar{\tau}_k^B)*G(\beta-\bar{\tau}_k^B)$$

```{code-cell}
#Polarisation P, bosonic

P_tau_b=G_tau_b*(basisf.u(beta-tausb.tau).T@G_l_f)

#(basisf.u(beta-tausb.tau).T@G_l_f) describes the evaluation on of G_l_f
#on sampling points (beta-tausb.tau)
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[1].plot(tausb.tau,P_tau_b.imag)
ax[1].set_title('Imaginary part of $P(tau)')
ax[0].plot(tausb.tau,P_tau_b.real)
ax[0].set_title('Real part of $P(tau)$')
```

The same way we transformed $\hat G^F(i\bar{\omega}_n^F)$ into $G^B(\bar\tau_k^B)$ we are now able to transform $P^B(\bar\tau_k^B)$ into $\hat P^B(i\bar{\omega}_n^B)$ again using least Square fitting and evaluating $P_l^B$ on the given sampling points (Matsubara frequencies).

+++

$$P_l^B$$

```{code-cell}
#P_l, bosonic

P_l_b=tausb.fit(P_tau_b)
```

```{code-cell}
pl.semilogy(np.abs(P_l_b),'+-')
pl.title('$P_l^B$ over l ')
pl.xticks(range(0,20))
pl.xlabel('l')
```

$$\hat{P}(i\bar{\omega}_k^B)= \sum_{l=0}^{N-1} P_l^BU_l^B(i\bar{\omega}_n^B)$$

```{code-cell}
#P_iw, bosonic

P_iw_b=np.einsum('lv,l->v', basisb.uhat(matsb.wn),P_l_b)
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[1].plot(matsb.wn,P_iw_b.imag)
ax[1].set_title('Imaginary part of $P(i\omega)$')
ax[0].plot(matsb.wn,P_iw_b.real)
ax[0].set_title('Real part of $P(i\omega)$')
```

Following $P$ we will now calculate the Screened Interaction $W$ with the following formula: 
 
$$\hat{{W}}(i\bar{\omega}_n^B)= U + U\hat{P}(i\bar{\omega}_n^B)\hat{{W}}(i\bar{\omega}_n^B)$$

or

$$\hat{{W}}(i\bar{\omega}_n^B)= \frac{U}{1-U\hat{P}(i\bar{\omega}_n^B)}$$

This equation has the exact same form as the Dyson equation but instead of connecting the Green's functions and the self energy it connects the Screened Coulomb Interaction $W$ and the polarization operator $P$. Because $U$ is the bare Coulomb interaction and $P$ is an object containing all irreducible processes (meaning in Feynmann diagrams cutting an interaction line U does not result in two sides) $W$ can be understood as the sum of all possible Feymann diagrams with an interaction U between its two parts.[5][11] This becomes quite intuitive if we look at 

$$W= U + UPW= U+UP(U+UP(U+UP(...)))=U+UPU+...$$

```{code-cell}
#W_iw, bosonic  

U=1/10

W_iw_b_U=U/(1-(U*P_iw_b))

W_iw_b=W_iw_b_U-U

#W_iw_b is the part depending on the frequency, any further calculations 
#will be done using this and not W_iw_b_U
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[1].plot(matsb.wn,W_iw_b.imag)
ax[1].set_title('Imaginary part of $W(i\omega)$')
ax[0].plot(matsb.wn,W_iw_b.real)
ax[0].set_title('Real part of $W(i\omega)$')
```

$$W_l^B$$

```{code-cell}
#W_l, bosonic

W_l_b=matsb.fit(W_iw_b)
```

```{code-cell}
pl.semilogy(np.abs(W_l_b),'+-')
pl.title('$W_l^B$ over l ')
pl.xticks(range(0,20))
pl.xlabel('l')
```

In the next step we are changing back into fermionic statistics: 

$${W}(\bar{\tau}_k^F)= \sum_{l=0}^{N-1} W_l ^BU_l^B(\bar{\tau}_k^F) $$

```{code-cell}
#W_tau_f, fermionic

W_tau_f=np.einsum('lt,l->t',basisb.u(tausf.tau),W_l_b)
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[1].plot(tausf.tau,W_tau_f.imag)
ax[1].set_title('Real part of $W(tau)$')
ax[0].plot(tausf.tau,W_tau_f.real)
ax[0].set_title('Imaginary part of $W(tau)$')
```

After changing back into fermionic statistics the next quantity that is being dealt with is the so-called self energy $\Sigma$. When we first introduced the Dyson equation in the first chapter, we saw $V$ as some sort of pertubation. The self-energy is very much alike as it describes the correlation effects of a many-body system. Here we will calculate $\tilde{\Sigma}(\bar{\tau}_k^F)$ using $\Sigma^{GW}$ so that

$${\Sigma}(\bar{\tau}_k^F)=-G(\bar{\tau}_k^F)*{W}(\bar{\tau}_k^F).$$

This can again be rewritten into its equivalent form

$$\Sigma=-G(U+UP(U+UP(...)))= -(GU+GUPU+GUPUPU+...)$$

```{code-cell}
#E_tau , fermionic   

E_tau_f=G_tau_f*W_tau_f
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[1].plot(tausf.tau,E_tau_f.imag)
ax[1].set_title('Imaginary part of $\Sigma(tau)$')
ax[0].plot(tausf.tau,E_tau_f.real)
ax[0].set_title('Real part of $\Sigma(tau)$')
```

$${\Sigma}_l^F$$

```{code-cell}
#E_l, fermionic

E_l_f=tausf.fit(E_tau_f)
```

```{code-cell}
pl.semilogy(np.abs(E_l_f),'+-')
pl.title('$\Sigma_l^F$ over l ')
pl.xticks(range(0,20))
pl.xlabel('l')
```

We will calculate

$$\hat{\Sigma}(i\bar{\omega}_n^F)=-UG(\beta)+\sum_{l=0}^{N-1}{\Sigma}_l^F \hat{U}_l^F(i\bar{\omega}_n^F)$$

with $UG(\beta)$ as the so called Hartee Fock term.

```{code-cell}
#E_iw, fermionic

E_iw_f=np.einsum('lw,l->w',basisf.uhat(matsf.wn),E_l_f)-U*(basisf.u(beta).T@G_l_f)
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[0].plot(matsf.wn,E_iw_f.imag)
ax[0].set_title('Imaginary part of $\Sigma(i\omega)$')
ax[1].plot(matsf.wn,E_iw_f.real)
ax[1].set_title('Real part of $\Sigma(i\omega)$')
```

Last but not least we will obtain $G^F(i\bar{\omega}_n^F)$ through the Dyson Equation of GW:

$$\hat{G}(i\bar{\omega}_k^F)=G_0(i\bar{\omega}_k^F)+G_0(i\bar{\omega}_k^F)\hat{\Sigma}(i\bar{\omega}_n^F)\hat{G}(i\bar{\omega}_k^F)$$

or

$$\hat{G}(i\bar{\omega}_k^F)=\frac{1}{G_0(i\bar{\omega}_k^F)^{-1}-\hat{\Sigma}(i\bar{\omega}_n^F)}$$

We have now calculated the full Green's function of the system with $\Sigma$ as a set of (one-particle) irreducible processes connected by $G_0$.[11]

```{code-cell}
#G_iw, fermonic     -> Dyson Equation

G_iw_f=((G_iw_0)**-1-E_iw_f)**-1        
```

```{code-cell}
fig, ax = pl.subplots(1,2, figsize=(12,5) )
ax[0].plot(matsf.wn,G_iw_f.imag)
ax[0].set_title('Imaginary part of $G(i\omega)')
ax[1].plot(matsf.wn,G_iw_f.real)
ax[1].set_title('Real part of $G(i\omega)')
```

```{code-cell}
print(G_iw_f)
```

With this we have completed one full iteration of GF2 & GW. Repeating this over and over again will show convergence and thus give the desired approximation to the self-energy.

+++

Sources
---------

+++

[1] H.Bruus, Karsten Flensberg,_Many-body Quantum Theory in Condensed    Matter  Physics_ (Oxford University Press Inc., New York, 2004)<br>
[2] J.W.Negele, H.Orland, _Quantum Many-Particle Systems_ (Westview Press, 1998)<br>
[3] G.D.Mahan,_Many-Particle Physics_ (Kluwer Academic/Plenum Publishers, New York, 2000)<br>
[4] P.C.Hansen,_Dicrete Inverse Problems_(Society for Industrial and Applied Mathematics, Philadelphia,2010)<br>
[5] J.Li, M.Wallerberger, N.Chikano, C.Yeh, E.Gull and H.Shinaoka, Sparse sampling approach to efficient ab initio calculations at finite temperature, Phys. Rev. B 101, 035144 (2020)<br>
[6] M.Wallerberger, H. Shinaoka, A.Kauch, Solving the Bethe–Salpeter equation with exponential convergence, Phys. Rev. Research 3, 033168 (2021)<br>
[7] H. Shinaoka, N. Chikano, E. Gull, J. Li, T. Nomoto, J. Otsuki, M. Wallerberger, T. Wang, K. Yoshimi, Efficient ab initio many-body calculations based on sparse modeling of Matsubara Green’s function (2021)<br>
[8] F. Aryasetiawan and O. Gunnarsson, The GW method, Rep. Prog. Phys. 61 (1998) 237<br>
[9] Matteo Gatti, _The GW approximation_, TDDFT school, Benasque (2014)<br>
[10]C.Friedrich, _Many-Body Pertubation Theory;The GW approximation_, Peter Grünberg Institut and Institute for Advanced Simulation, Jülich <br>
[11] K.. Held, C. Taranto, G. Rohringer, and A. Toschi, _Hedin Equations, GW, GW+DMFT, and All That_ (Modeling and Simulation Vol. 1, Forschungszentrum Jülich, 2011 )

```{code-cell}

```
