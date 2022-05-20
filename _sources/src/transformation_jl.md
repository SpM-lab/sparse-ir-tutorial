---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Julia 1.7
  language: julia
  name: julia-1.7
---

+++ {"tags": []}

# Transformation from/to IR

+++

## Poles

```{code-cell}
using SparseIR
using Plots
gr() # USE GR backend
using OMEinsum
using LaTeXStrings

#import PyPlot #To be removed


beta = 15.0
wmax = 10.0
basis_b = FiniteTempBasis(boson, beta, wmax, 1e-7)

coeff = [1.0]
omega_p = [0.1]

rhol_pole = ein"lp,p->l"(
    basis_b.v(omega_p),
    coeff ./ tanh.(0.5*beta*omega_p)
)
gl_pole = - basis_b.s .* rhol_pole

plot([abs.(rhol_pole), abs.(gl_pole)], yaxis=:log, label=[latexstring("|\\rho_l|") latexstring("|g_l|")], xlabel=latexstring("l"), ylims=(1e-5, 1e+1))
```

```{code-cell}
sp = SparsePoleRepresentation(basis_b, omega_p)
gl_pole2 = to_IR(sp, coeff)

plot(
    [abs.(gl_pole2), abs.(gl_pole)],
    label=["from SPR" "exact"],
    yaxis=:log,
    xlabel=latexstring("l"), ylabel=latexstring("g_l"), ylims=(1e-5, 1e+1), marker=[:cross :circle])
```

## From smooth spectral function

```{code-cell}
# Three Gaussian peaks (normalized to 1)
gaussian(x, mu, sigma) = exp(-((x-mu)/sigma)^2)/(sqrt(π)*sigma)

rho(omega) = 0.2*gaussian(omega, 0.0, 0.15) + 
    0.4*gaussian(omega, 1.0, 0.8) + 0.4*gaussian(omega, -1.0, 0.8)

omegas = LinRange(-5, 5, 1000)
plot(omegas, rho.(omegas), xlabel=latexstring("\\omega"), ylabel=latexstring("\\rho(\\omega)"), label="")
```

```{code-cell}
beta = 10.0
wmax = 10.0
basis = FiniteTempBasis(fermion, beta, wmax, 1e-7)

rhol = [overlap(basis.v[l], rho) for l in 1:length(basis)]
gl = - basis.s .* rhol

plot([abs.(rhol), abs.(gl)], yaxis=:log, ylims=(1e-5,1), marker=[:circle :diamond], line=nothing, xlabel=latexstring("l"), label=[latexstring("\\rho_l") latexstring("|g_l|")])
```

```{code-cell}
rho_omega_reconst = transpose(basis.v(omegas)) * rhol

plot(omegas, rho_omega_reconst, xlabel=latexstring("\\omega"), ylabel=latexstring("\\rho(\\omega)"), label="")
```

## From IR to imaginary time

```{code-cell}
taus = collect(LinRange(0, beta, 1000))
gtau1 = transpose(basis.u(taus)) * gl

plot(taus, gtau1, xlabel=latexstring("\\tau"), ylabel=latexstring("G(\\tau)"), label=nothing)
```

```{code-cell}
smpl = TauSampling(basis, taus)
gtau2 = evaluate(smpl, gl)
plot(taus, gtau2, xlabel=latexstring("\\tau"), ylabel=latexstring("G(\\tau)"), label=nothing, marker=:cross)
```

## From full imaginary-time data

```{code-cell}
function eval_gtau(taus)
    uval = basis.u(taus) #(nl, ntau)
    return transpose(uval) * gl
end

gl_reconst = [overlap(basis.u[l], eval_gtau) for l in 1:length(basis)]

ls = collect(0:length(basis)-1)
plot(
    ls[1:2:end],
    [abs.(gl_reconst[1:2:end]), abs.(gl[1:2:end]), abs.(gl_reconst - gl)[1:2:end]], xlabel=latexstring("l"), label=["reconstructed" "exact" "error"], marker=[:+ :x :none], markersize=10, yaxis=:log)
```

```{code-cell}

```
