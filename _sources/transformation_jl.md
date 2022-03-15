---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.7
kernelspec:
  display_name: Julia 1.7.2
  language: julia
  name: julia-1.7
---

+++ {"tags": []}

# Transformation from/to IR

+++

## Poles

```{code-cell}
using SparseIR
#ENV["MPLBACKEND"]="tkagg"
using PyPlot
using OMEinsum


beta = 15
wmax = 10
basis_b = FiniteTempBasis(boson, beta, wmax, 1e-10)

coeff = [1.0]
omega_p = [0.1]

rhol_pole = ein"lp,p->l"(
    basis_b.v(omega_p),
    coeff ./ tanh.(0.5*beta*omega_p)
)
gl_pole = - basis_b.s .* rhol_pole

plt.semilogy(abs.(rhol_pole), marker="o", label=L"|\rho_l|")
plt.semilogy(abs.(gl_pole), marker="x", label=L"|g_l|")

plt.xlabel(L"$l$")
plt.ylim([1e-5, 1e+1])
plt.legend()
;
```

```{code-cell}
sp = SparsePoleRepresentation(basis_b, omega_p)
gl_pole2 = to_IR(sp, coeff)

plt.semilogy(abs.(gl_pole2), marker="x", label=L"$|g_l|$ from SPR")
plt.semilogy(abs.(gl_pole), marker="x", label=L"$|g_l|$")

plt.xlabel(L"$l$")
plt.ylim([1e-5, 1e+1])
plt.legend()
#plt.show()
```

## From smooth spectral function

```{code-cell}
# Three Gaussian peaks (normalized to 1)
gaussian(x, mu, sigma) = exp(-((x-mu)/sigma)^2)/(sqrt(π)*sigma)

rho(omega) = 0.2*gaussian(omega, 0.0, 0.15) + 
    0.4*gaussian(omega, 1.0, 0.8) + 0.4*gaussian(omega, -1.0, 0.8)

omegas = LinRange(-5, 5, 1000)
plt.xlabel(L"\omega")
plt.ylabel(L"\rho(\omega)")
plt.plot(omegas, rho.(omegas))
#plt.show()
```

```{code-cell}
beta = 10
wmax = 10
basis = FiniteTempBasis(fermion, beta, wmax, 1e-10)

rhol = overlap(basis.v, rho)
gl = - basis.s .* rhol

plt.semilogy(abs.(rhol), marker="o", label=L"|\rho_l|")
plt.semilogy(abs.(gl), marker="x", label=L"|g_l|")
plt.xlabel(L"l")
plt.ylim([1e-5, 1])
plt.legend()
#plt.show()
```

```{code-cell}
rho_omgea_reconst = transpose(basis.v(omegas)) * rhol

plt.xlabel(L"\omega")
plt.ylabel(L"\rho(\omega)")
plt.plot(omegas, rho.(omegas))
#plt.show()
```

## From IR to imaginary time

```{code-cell}
taus = collect(LinRange(0, beta, 1000))
gtau1 = transpose(basis.u(taus)) * gl
plt.plot(taus, gtau1)
plt.xlabel(L"\tau")
plt.ylabel(L"G(\tau)")
#plt.show()
```

```{code-cell}
smpl = TauSampling(basis, taus)
gtau2 = evaluate(smpl, gl)
plt.plot(taus, gtau1)
plt.xlabel(L"\tau")
plt.ylabel(L"G(\tau)")
#plt.show()
```

## From full imaginary-time data

```{code-cell}
function eval_gtau(taus)
    uval = basis.u(taus) #(nl, ntau)
    return transpose(uval) * gl
end

gl_reconst = overlap(basis.u, eval_gtau)

plt.semilogy(abs.(gl_reconst), label="reconstructed", marker="o")
plt.semilogy(abs.(gl), label="exact", marker="x")
plt.semilogy(abs.(gl_reconst - gl), label="error", marker="")
plt.xlabel(L"l")
plt.ylabel(L"|g_l|")
plt.ylim([1e-20, 1])
plt.legend()
#plt.show()
```

```{code-cell}

```
