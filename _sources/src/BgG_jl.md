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

# Bogoliubov-de Gennes equations in real space

Author: [Yuki Nagai](mailto:nagai.yuki@jaea.go.jp)

**This notebook is expensive. Running this notebook takes a few minutes.**

We consider a superconductivity in a mean-field level.

## Basic theory of superconductivity in real space
Let us consider a Hamiltonian for a fermion system given as $H = \psi^{\dagger} \hat{\cal H} \psi/2$. 
The column vector $\psi$ is composed of $N$ fermionic annihilation $c_i$ and creation operators $c_i^{\dagger}$ $(i = 1,2,\cdots,N)$, $\psi = (\{ c_i \}^T, \{ c_i^{\dagger} \}^T)$, where $\{ c_i \} = (c_1,c_2,\cdots,c_N)^T$ and $\{ c_i^{\dagger} \} = (c_1^{\dagger},c_2^{\dagger},\cdots,c_N^{\dagger})$. The row vector $\psi^{\dagger}$ is also defined as $\psi^{\dagger} = (\{c_i^{\dagger} \}^T, \{ c_i \}^T)$.
The subscription $i$ in $c_i$ or $c_i^{\dagger}$ indicates a quantum index depending on spatial site, spin, orbital, etc. The "Hamiltonian" matrix $\hat{\cal H}$ is a $2N \times 2N$ Hermitian matrix given as 

$$
\hat{\cal H} = \left( \begin{matrix} \hat{H} & \hat{\Delta} \\
\hat{\Delta}^{\dagger} & - \hat{H}^{\ast}
\end{matrix}
\right),
$$

where $\hat{H}$ is a normal state Hamiltonian and $\hat{\Delta}$ is a superconducting order parameter.

+++

The BdG equations are regarded as the eigenvalue equation with respect to $\hat{\cal H}$ expressed as 

$$
\hat{\cal H} \vec{f}_{\gamma} = \epsilon_{\gamma} \vec{f}_{\gamma}
$$

$$
 \vec{f}_{\gamma} = \left(  \begin{matrix}
 \vec{u}_{\gamma} \\
 \vec{v}_{\gamma}
 \end{matrix}
 \right)
$$

The column vectors $\vec{u}_{\gamma}$ and $\vec{v}_{\gamma}$ are $N$-component complex vetors. 
To solve the BdG equations is equivalent to diagonalization of $\hat{\cal H}$ with a unitary matrix $\hat{U}$, 

$$
\hat{U}^{\dagger} \hat{\cal H} \hat{U} = \hat{D}, \: \: \: \hat{D} = {\rm diag} (\epsilon_1,\epsilon_2,\cdots,\epsilon_{2N})
$$

+++

Mean-fields are calculated by one-particle Green's functions. 
In the mean-field framework, the $2N \times 2N$ matrix Green's function with a complex frequency is defined as 

$$
\hat{G}(z) = [z \hat{I} - \hat{\cal H}]^{-1}.
$$

With the unitary matrix $\hat{U}$, each component of $\hat{G}(z)$ is epressed as 

$$
G_{\alpha \beta}(z) = \sum_{\gamma=1}^{2N} U_{\alpha \gamma} U_{\beta \gamma}^{\ast} \frac{1}{z - \epsilon_{\gamma}}
$$

If we set $z = i \omega_n$ with the Matsubara frequency $\omega_n = (2n+1)\pi T$, the above formula corresponds to Matsubara temperature Green's function. 
The regarded and advanced Green's functions are, respectively, defined as 

$$
\hat{G}^{\rm R}(\omega) = \lim_{\tau \rightarrow 0+} \hat{G}(\omega + i \eta)
$$

$$
\hat{G}^{\rm A}(\omega) = \lim_{\tau \rightarrow 0+} \hat{G}(\omega - i \eta)
$$

In order to obtain physical quantities (e.g., density of states) from Green's functions, we introduce the following useful $2N$-component unit-vectors $\vec{e}(i)$ and $\vec{h}(i)$ $(1 \le i \le N)$, which are, respectively, defined as

$$
[\vec{e}(i)]_{\gamma} = \delta_{i,\gamma}, \: [\vec{h}(i)]_{\gamma} = \delta_{i+N,\gamma}
$$

For example, the local density of states with respect to the site $i$ is given as 

$$
N(\omega,i) = -\frac{1}{2\pi i} \vec{e}(i)^T \hat{d}(\omega) \vec{e}(i),
$$

with the use of the spectral function $\hat{d}(\omega)$ defined as 

$$
[\hat{d}(\omega)]_{\alpha \beta} \equiv [\hat{G}^{\rm R}(\omega) - \hat{G}^{\rm A}(\omega)]_{\alpha \beta} = - 2\pi i \sum_{\gamma=1}^{2N} U_{\alpha \gamma} U_{\beta \gamma}^{\ast} \delta(\omega - \epsilon_{\gamma})
$$

Two types of mean-fields $\langle c_i^{\dagger} c_j \rangle$ and $\langle c_i c_j \rangle$ can be expressed as 

$$
\langle c_i^{\dagger} c_j \rangle = - \frac{1}{2 \pi i} \int_{-\infty}^{\infty} d\omega f(\omega) \vec{e}(j)^T \hat{d}(\omega) \vec{e}(i), 
$$

$$
\langle c_i c_j \rangle = - \frac{1}{2 \pi i} \int_{-\infty}^{\infty} d\omega f(\omega) \vec{e}(j)^T \hat{d}(\omega) \vec{h}(i), 
$$

with $f(x) \equiv 1/(e^{x/T} + 1)$. 

+++

## Matsubara frequency representation
With the use of the analytic continuation, the mean-fields are rewritten as 

$$
\langle c_i^{\dagger} c_j \rangle = T \sum_{n=-\infty}^{\infty} \vec{e}(j)^T \hat{G}(i \omega_n) \vec{e}(i)
$$

$$
\langle c_i c_j \rangle = T \sum_{n=-\infty}^{\infty} \vec{e}(j)^T \hat{G}(i \omega_n) \vec{h}(i)
$$

By solving the linear equations defined as 

$$
(i \omega_n \hat{I} - \hat{\cal H}) \vec{x}(i,\omega_n) = \vec{h}(i), 
$$

the superconducting mean-field is expressed as 

$$
\langle c_i c_j \rangle = T \sum_{n = - \infty}^{\infty} \hat{e}(j)^T \vec{x}(i,\omega_n).
$$

For example, s-wave superconducting order parameter is defined as 

$$
\Delta_i = U \langle c_i c_i \rangle,
$$

The self-consistent cycle has the schematic form 

$$ 
\hat{\Delta} \rightarrow \hat{G}(i \omega_n) \rightarrow \langle c_i c_i \rangle \rightarrow \hat{\Delta}
$$



+++

## Reduced-shifted Conjugation Gradient method
The reduced-shifted conugate-gradient (RSCG) method is numerically efficient to calculate a matrix element of a Green's function defined as a resolvent of a Hamiltonian operator, by solving linear equations with desired accuracy. The matrix elements with different frequencies are simultaneously obtained. The details are described in J. Phys. Soc. Jpn. `86`, 014708 (2017). The RSCG method is provided in RSCG.jl package. 

+++

## IR basis
The superconducting mean-fielads are expressed as 

$$
\langle c_i c_j \rangle =  \vec{e}(j)^T \hat{G}(\tau = \beta) \vec{h}(i) = \sum_{l=0}^{N_{\rm IR}-1} U_l(\beta) \vec{e}(j)^T \hat{G}_l \vec{h}(i)
$$

Since the Matsubara Green's function is expressed as 

$$
G(i \omega_n) = \sum_{l=0}^{N_{\rm IR}-1} G_l U_l(i \omega_n),
$$

one can obtain $G_l$ by fitting the above equation. 

+++

# Two-dimensional s-wave superconductor on the square lattice

The normal state Hamiltonian on the tight-binding model is expressed as 

$$
\sum_{i,j} t_{ij} c_i^{\dagger} c_j - \mu \sum_i c_i^{\dagger} c_i
$$
Here, $t_{ij} = -t$ if $j$ is a nearest neighbor for $i$. 
The matrix form of the above Hamiltonain is expressed as 

$$
\hat{H} = -\mu \hat{I} -t (\hat{T}_{+x} + \hat{T}_{-x} + \hat{T}_{+y} + \hat{T}_{-y}),
$$

where $\hat{T}_{+x}$ is a hopping matrix with respect to $+x$ direction. 


+++

In Julia code, this Hamiltonian matrix is defined as 

```{code-cell}

```

```{code-cell}
# Check if a given function called with given types is type stable
function typestable(@nospecialize(f), @nospecialize(t))
    v = code_typed(f, t)
    stable = true
    for vi in v
        for (name, ty) in zip(vi[1].slotnames, vi[1].slottypes)
            !(ty isa Type) && continue
            if ty === Any
                stable = false
                println("Type instability is detected! the variable is $(name) ::$ty")
            end
        end
    end
    return stable
end

using SparseArrays
using LinearAlgebra
function make_x_plushop(Nx,Ny,BC)
    N = Nx*Ny
    Txhop = spzeros(Int64,N,N)
    for ix=1:Nx
        for iy=1:Ny
            i = (iy-1)*Nx + ix
            jx = ix + 1
            jy = iy
            if BC == "PBC"
                jx += ifelse(jx > Nx,-Nx,0)
            elseif BC == "OBC"
            else
                error("BC = $BC is not supported")
            end
            if 1 <= jx <= Nx
                j = (jy-1)*Nx + jx
                Txhop[i,j] = 1
            end
        end
    end
    return Txhop
end

function make_x_minushop(Nx,Ny,BC)
    N = Nx*Ny
    Txhop = spzeros(Int64,N,N)
    for ix=1:Nx
        for iy=1:Ny
            i = (iy-1)*Nx + ix
            jx = ix - 1
            jy = iy
            if BC == "PBC"
                jx += ifelse(jx < 1,Nx,0)
            elseif BC == "OBC"
            else
                error("BC = $BC is not supported")
            end
            if 1 <= jx <= Nx
                j = (jy-1)*Nx + jx
                Txhop[i,j] = 1
            end
        end
    end
    return Txhop
end

function make_y_plushop(Nx,Ny,BC)
    N = Nx*Ny
    Tyhop = spzeros(Int64,N,N)
    for ix=1:Nx
        for iy=1:Ny
            i = (iy-1)*Nx + ix
            jx = ix 
            jy = iy + 1
            if BC == "PBC"
                jy += ifelse(jy > Ny,-Ny,0)
            elseif BC == "OBC"
            else
                error("BC = $BC is not supported")
            end
            if 1 <= jy <= Ny
                j = (jy-1)*Nx + jx
                Tyhop[i,j] = 1
            end
        end
    end
    return Tyhop
end

function make_y_minushop(Nx,Ny,BC)
    N = Nx*Ny
    Tyhop = spzeros(Int64,N,N)
    for ix=1:Nx
        for iy=1:Ny
            i = (iy-1)*Nx + ix
            jx = ix 
            jy = iy - 1
            if BC == "PBC"
                jy += ifelse(jy < 1,Ny,0)
            elseif BC == "OBC"
            else
                error("BC = $BC is not supported")
            end
            if 1 <= jy <= Ny
                j = (jy-1)*Nx + jx
                Tyhop[i,j] = 1
            end
        end
    end
    return Tyhop
end

function make_H_normal(Nx,Ny,μ,BC)
    N = Nx*Ny
    Tx_plushop = make_x_plushop(Nx,Ny,BC)
    Tx_minushop = make_x_minushop(Nx,Ny,BC)
    Ty_plushop = make_y_plushop(Nx,Ny,BC)
    Ty_minushop = make_y_minushop(Nx,Ny,BC)
    HN = sparse(I,N,N)*(-μ)
    t = 1.0
    
    HN += -t*(Tx_plushop + Tx_minushop + Ty_plushop + Ty_minushop)
    return HN
end

@assert typestable(make_x_plushop,(Int64,Int64,String))
@assert typestable(make_x_minushop,(Int64,Int64,String))
@assert typestable(make_y_plushop,(Int64,Int64,String))
@assert typestable(make_y_minushop,(Int64,Int64,String))
```

The onsite superconducting s-wave order parameter is defined as 

```{code-cell}
function make_Δ(Δ)
    Nx,Ny = size(Δ)
    N = Nx*Ny
    Δmat = spzeros(ComplexF64,N,N)
    for ix=1:Nx
        for iy=1:Ny
            i = (iy-1)*Nx + ix
            Δmat[i,i] = Δ[ix,iy]
        end
    end
    return Δmat
end
@assert typestable(make_Δ,(Matrix{ComplexF64},))
```

The BdG Hamiltonian is defined as 

```{code-cell}
function make_H_sc(Nx,Ny,μ,Δ,BC)
    HN = make_H_normal(Nx,Ny,μ,BC)
    matΔ = make_Δ(Δ)
    H = [
        HN matΔ
        matΔ' -conj.(HN)
    ]
    return H
end
@assert typestable(make_H_sc,(Int64,Int64,Float64,Matrix{ComplexF64},String))
```

We have to update the order parameters in the self-consistent cycle.
The update method is implemented as 

```{code-cell}
function update_H_sc!(H,Δ)
    matΔ = make_Δ(Δ)
    Nx,Ny = size(Δ)
    N = Nx*Ny
    H[1:N,1+N:2N] = matΔ
    H[1+N:2N,1:N] = matΔ'
end

@assert typestable(update_H_sc!,(AbstractMatrix{ComplexF64},Matrix{ComplexF64}))
```

## Conventional method with Matsubara frequencies
If one wants to use the Matsubara frequencies, one has to introduce the cutoff of the frequencies. Then the mean-fields are expressed as 

$$
\langle c_i c_j \rangle = T \sum_{n = - n_c}^{n_c} \hat{e}(j)^T \vec{x}(i,\omega_n).
$$

Here, $n_c$ is a cutoff value. 
With the use of the maximum frequency $\omega_c$, we calculate the set of Matsubara frequencies defined as 

```{code-cell}
function calc_ωn(T,ωc)
    M = Int((round(ωc/(T*π)))/2-1)
    println("num. of Matsubara freq: ",2M)
    ωn = zeros(ComplexF64,2M)
    for n=1:2M
        ωn[n] = π*T*(2.0*(n-M-1)+1)*im
    end
    return ωn
end
@assert typestable(calc_ωn,(Float64,Float64))
```

```{code-cell}
ωc = 100π
T = 0.01
ωn =  calc_ωn(T,ωc)
```

With the use of RSCG.jl, the mean-fields are calculated by 

```{code-cell}
using Pkg
Pkg.add(PackageSpec(name="RSCG",rev = "master"))
using RSCG
```

```{code-cell}
function calc_Δi!(i,N,H,Δold,T,U,ωn;mixratio = 0.5)
    j = i + N
    Gij = greensfunctions(i,j,ωn,H) 
    Δi = U*T*sum(Gij)
    Δi = (1-mixratio)*Δold[i] + mixratio*Δi   
    return Δi
end

@assert typestable(calc_Δi!,(Int64,Int64,AbstractMatrix{ComplexF64},Matrix{ComplexF64},Float64,Float64,Vector{ComplexF64}))

function calc_Δ!(Δnew,H,Δold,T,U,ωn;mixratio = 0.5)
    Nx,Ny = size(Δold)
    N = Nx*Ny
    map!(i -> calc_Δi!(i,N,H,Δold,T,U,ωn,mixratio = mixratio),Δnew,1:N) #If you use pmap! instead of map!, you can do the parallel computation.
    return 
end

@assert typestable(calc_Δi!,(AbstractMatrix{ComplexF64},Matrix{ComplexF64},Float64,Float64,Vector{ComplexF64}))
```

Then, we can construct the BdG Hamiltonian matrix. 

```{code-cell}
Nx = 10
Ny = 10
Δ = ones(ComplexF64,Nx,Ny)
Δold = copy(Δ)
Δnew = zero(Δ)
BC = "OBC" #open boundary condition
#BC = "PBC" #periodic boundary condition
U  =-2
μ = -0.2

H = make_H_sc(Nx,Ny,μ,Δ,BC)
```

The self-consistent loop is defined as 

```{code-cell}
itemax = 100
ix = Nx ÷ 2
iy = Ny ÷ 2
for ite = 1:itemax
    calc_Δ!(Δnew,H,Δold,T,U,ωn)
    update_H_sc!(H,Δnew)
    eps = sum(abs.(Δnew-Δold))/sum(abs.(Δold))
    println("$ite $eps ",Δnew[ix,iy])
    Δold .= Δnew
    if eps < 1e-4
        break
    end
end
```

## IR based method
IR based method for BdG equations is shown in arXiv:2111.13288. In this paper, we use the ultra-fast method (J. Phys. Soc. Jpn. 89, 074703 (2020)). In this notebook, we do not use the ultra-fast method for simplicity. 

We use the SparseIR package. The sampling Matsubara frequencies are obtained as 

```{code-cell}
using SparseIR
wmax = 10.0
beta = 1/T
basis = FiniteTempBasis(fermion, beta, wmax, 1e-7)
smpl = MatsubaraSampling(basis)
ωn_s = smpl.sampling_points .* (im * π/beta)
println("num. of Matsubara freqs. ", length(ωn_s))
smpl_beta = TauSampling(basis, [beta])

```

The mean-fields are obtained by 

```{code-cell}
function fit_ir(Gij,smpl_Matsubara,smpl_beta)
    gl = fit(smpl_Matsubara, Gij)
    G0 = evaluate(smpl_beta, gl)
    return -G0[1]
end

@assert typestable(fit_ir,(Vector{ComplexF64},typeof(smpl),typeof(smpl_beta)))

function calc_Δi_ir!(i,N,H,Δold,T,U,ωn,smpl_Matsubara,smpl_beta;mixratio = 0.5)
    j = i + N
    Gij = greensfunctions(i,j,ωn,H) 
    G0 = fit_ir(Gij,smpl_Matsubara,smpl_beta)            
    Δi = U*G0
    Δi = (1-mixratio)*Δold[i] + mixratio*Δi   
    return Δi
end

@assert typestable(calc_Δi_ir!,(Int64,Int64,AbstractMatrix{ComplexF64},Matrix{ComplexF64},Float64,
            Float64,Vector{ComplexF64},typeof(smpl),typeof(smpl_beta)))


function calc_Δ_ir!(Δnew,H,Δold,T,U,ωn,smpl_Matsubara,smpl_beta;mixratio = 0.5)
    Nx,Ny = size(Δold)
    N = Nx*Ny
    map!(i -> calc_Δi_ir!(i,N,H,Δold,T,U,ωn,smpl_Matsubara,smpl_beta,mixratio = mixratio),Δnew,1:N) #If you use pmap! instead of map!, you can do the parallel computation.
    return
end

@assert typestable(calc_Δ_ir!,(Matrix{ComplexF64},AbstractMatrix{ComplexF64},Matrix{ComplexF64},Float64,
            Float64,Vector{ComplexF64},typeof(smpl),typeof(smpl_beta)))
```

Then, the self-consistent cycle is expressed as 

```{code-cell}
Nx = 10
Ny = 10
Δ = ones(ComplexF64,Nx,Ny)
Δold = copy(Δ)
Δnew = zero(Δ)
BC = "OBC" #open boundary condition
#BC = "PBC" #periodic boundary condition
U  =-2
μ = -0.2

H = make_H_sc(Nx,Ny,μ,Δ,BC)

for ite = 1:itemax
    calc_Δ_ir!(Δnew,H,Δold,T,U,ωn_s,smpl,smpl_beta)
    update_H_sc!(H,Δnew)
    
    eps = sum(abs.(Δnew-Δold))/sum(abs.(Δold))
    println("$ite $eps ",Δnew[ix,iy])
    Δold .= Δnew
    if eps < 1e-4
        break
    end
end
```

Let us show the superconducting order parameter in real space. Since we consider the open boundary condition, the mean-field becomes large at four corner due to an interference effect. 

```{code-cell}
using Plots
plot(1:Nx,1:Ny,abs.(Δnew),st=:surface,
    xlabel = "ix",
    ylabel = "iy",
    zlabel = "|Delta|")
```

We plot the local density of states at the center.

```{code-cell}
using Plots
M = 1000
σ = zeros(ComplexF64,M)
η = 0.05
σmin = -4.0 + im*η
σmax = 4.0+ im*η
for i=1:M
    σ[i] = (i-1)*(σmax-σmin)/(M-1) + σmin
end

i = (iy-1)*Nx + ix
j = i

Gij1 = greensfunctions(i,j,σ,H) 
plot(real.(σ),(-1/π)*imag.(Gij1),
    xlabel = "Energy [t]",
    ylabel = "Local DOS",)

```

# Two-dimensional topological s-wave superconductor on the square lattice
We consider topological superconductivity. 
The s-wave superconductor becomes topologically non-trivial when there are Zeeman magnetic fields and Rashba spin-orbit coupling. 
The Hamiltonian in momentum space is given as 

$$
\hat{\cal H}(\vec{k}) = 
\left(
    \begin{matrix}
    h_0(\vec{k}) & i \Delta \sigma_y \\
    - i \Delta^{\dagger} \sigma_y & - h_0^{\ast}(\vec{k})
\end{matrix}
\right), 
$$

with $ h_0(\vec{k}) = \epsilon(\vec{k}) - h \sigma_z + \alpha {\cal L}(\vec{k})$. 
The band dispersion is $\epsilon(\vec{k}) = -2 t (\cos k_x + \cos k_y) - \mu$. 
The Zeeman magnetic field is $h$. 
The Rashba spin-orbit term is described by $\alpha {\cal L}(\vec{k}) = \alpha (\sigma_x \sin k_y - \sigma_y \sin k_x)$. 
Here, $\sigma_{x,y,z}$ are the Pauli matrices. 

By Fourier transforming this Hamiltonian, we have the Hamiltonian defined in real space.

+++

We difine the Pauli matrices. 

```{code-cell}
const σ0 = [1 0
0 1]
const σx = [0 1
1 0]
const σy = [0 -im
im 0]
const σz = [1 0
0 -1]
```

We can construct the Hamiltonian with the use of the ```kron```: 

```{code-cell}
function make_Htsc_normal(Nx,Ny,μ,BC,h,α)
    N = Nx*Ny
    Tx_plushop = make_x_plushop(Nx,Ny,BC)
    Tx_minushop = make_x_minushop(Nx,Ny,BC)
    Ty_plushop = make_y_plushop(Nx,Ny,BC)
    Ty_minushop = make_y_minushop(Nx,Ny,BC)
    HN = kron(sparse(I,N,N)*(-μ),σ0) 
    HN += kron(sparse(I,N,N)*(-h),σz) #Zeeman magnetic field
    t = 1.0
    
    HN += kron(-t*(Tx_plushop + Tx_minushop + Ty_plushop + Ty_minushop),σ0)
    
    Hax = kron((α/(2im))*(Tx_plushop - Tx_minushop ) ,σy)
    HN += Hax 
    Hay = kron((α/(2im))*(Ty_plushop - Ty_minushop ) ,σx)
    HN += Hay 
    
    return HN
end

@assert typestable(make_Htsc_normal,(Int64,Int64,Float64,String,Float64,Float64))
```

We also define the superconducting order parameter: 

```{code-cell}

function make_Δtsc(Δ)
    Nx,Ny = size(Δ)
    N = Nx*Ny
    Δmat = spzeros(ComplexF64,N,N)
    for ix=1:Nx
        for iy=1:Ny
            i = (iy-1)*Nx + ix
            Δmat[i,i] = Δ[ix,iy]
        end
    end
    return kron(Δmat,im*σy)
end

@assert typestable(make_Δtsc,(Matrix{ComplexF64},))
```

Then, the superconducting BdG Hamiltonian is defined as 

```{code-cell}
function make_Htsc_sc(Nx,Ny,μ,Δ,BC,h,α)
    HN = make_Htsc_normal(Nx,Ny,μ,BC,h,α)
    matΔ = make_Δtsc(Δ)
    H = [
        HN matΔ
        matΔ' -conj.(HN)
    ]
    return H
end

@assert typestable( make_Htsc_sc,(Int64,Int64,Float64,Matrix{ComplexF64},String,Float64,Float64))
```

The update method is implimanted as 

```{code-cell}
function update_Htsc_sc!(H,Δ)
    matΔ = make_Δtsc(Δ)
    N,_ = size(matΔ)
    H[1:N,1+N:2N] = matΔ
    H[1+N:2N,1:N] = matΔ'
end

@assert typestable(update_Htsc_sc!,(AbstractMatrix{ComplexF64},Matrix{ComplexF64}))
```

The mean-fields are calculated as 

```{code-cell}
function calc_Δitsc_ir!(i,N,H,Δold,T,U,ωn,smpl_Matsubara,smpl_beta;mixratio = 0.5)
    ispin = 1
    ii = (i-1)*2 + ispin
    jspin = 2
    jj = (i-1)*2 + jspin + N
    
    Gij = greensfunctions(ii,jj,ωn,H) 
    G0 = fit_ir(Gij,smpl_Matsubara,smpl_beta)            
    Δi = U*G0
    Δi = (1-mixratio)*Δold[i] + mixratio*Δi   
    return Δi
end

@assert typestable(calc_Δitsc_ir!,(Int64,Int64,AbstractMatrix{ComplexF64},Matrix{ComplexF64},Float64,
            Float64,Vector{ComplexF64},typeof(smpl),typeof(smpl_beta)))


function calc_Δtsc_ir!(Δnew,H,Δold,T,U,ωn,smpl_Matsubara,smpl_beta,;mixratio = 0.5)
    Nx,Ny = size(Δold)
    N = Nx*Ny*2
    map!(i -> calc_Δitsc_ir!(i,N,H,Δold,T,U,ωn,smpl_Matsubara,smpl_beta,mixratio = mixratio),Δnew,1:Nx*Ny) #If you use pmap! instead of map!, you can do the parallel computation.
end

@assert typestable(calc_Δtsc_ir!,(Matrix{ComplexF64},AbstractMatrix{ComplexF64},Matrix{ComplexF64},Float64,
            Float64,Vector{ComplexF64},typeof(smpl),typeof(smpl_beta)))
```

Then, we can solve the gap equations, self-consistently. 

```{code-cell}
T = 0.01

beta = 1/T
wmax = 10.0

basis = FiniteTempBasis(fermion, beta, wmax, 1e-7)
smpl = MatsubaraSampling(basis)
ωn = smpl.sampling_points .* (im * π/beta)
println("num. of Matsubara freqs. ", length(ωn))
smpl_beta = TauSampling(basis, [beta])

U  =-5.6
itemax = 1000
μ = 3.5

Nx = 16
Ny = 16
Δ = 3*ones(ComplexF64,Nx,Ny)
Δold = copy(Δ)
Δnew = zero(Δ)
BC = "OBC"
h = 1
α = 1
Htsc =  make_Htsc_sc(Nx,Ny,μ,Δold,BC,h,α)

ix = Nx ÷ 2
iy = Ny ÷ 2
i = (iy-1)*Nx + ix
j = i


for ite = 1:itemax
    calc_Δtsc_ir!(Δnew,Htsc,Δold,T,U,ωn,smpl,smpl_beta,mixratio =1)
    update_Htsc_sc!(Htsc,Δnew)
    
    eps = sum(abs.(Δnew-Δold))/sum(abs.(Δold))
    println("$ite $eps ",Δnew[ix,iy]," ave: ",sum(abs.(Δnew))/length(Δnew))
    
    Δold .= Δnew
    if eps < 1e-3
        break
    end
end
```

The superconducting order parameter is shown as 

```{code-cell}
using Plots
plot(1:Nx,1:Ny,abs.(Δnew),st=:surface,
    xlabel = "ix",
    ylabel = "iy",
    zlabel = "|Delta|")
```

We can see that the order parametere oscillates due to the Fridel oscillation. 

+++

The local density of states at the center is calculated as 

```{code-cell}
using Plots
M = 1000
σ = zeros(ComplexF64,M)
η = 0.01
σmin = -4.0 + im*η
σmax = 4.0+ im*η
for i=1:M
    σ[i] = (i-1)*(σmax-σmin)/(M-1) + σmin
end


ix = Nx ÷ 2
iy = Ny ÷ 2
ispin  =1

i = (( iy-1)*Nx + ix-1)*2 + ispin
j = i
Gij1 = greensfunctions(i,j,σ,Htsc) 

ispin  =2
i = (( iy-1)*Nx + ix-1)*2 + ispin
j = i
Gij2 = greensfunctions(i,j,σ,Htsc) 
plot(real.(σ),(-1/π)*imag.(Gij1 .+ Gij2),
    xlabel = "Energy [t]",
    ylabel = "Local DOS",)

```

At the edge, there are mid-gap states since this is topologically non-trivial state. We can see the local density of states at the edge. 

```{code-cell}

ix = 1
iy = Ny ÷ 2
ispin  =1

i = (( iy-1)*Nx + ix-1)*2 + ispin
j = i
Gij1 = greensfunctions(i,j,σ,Htsc) 

ispin  =2
i = (( iy-1)*Nx + ix-1)*2 + ispin
j = i
Gij2 = greensfunctions(i,j,σ,Htsc) 
plot(real.(σ),(-1/π)*imag.(Gij1 .+ Gij2),
    xlabel = "Energy [t]",
    ylabel = "Local DOS",)
```

We can see there are mid-gap states.

```{code-cell}

```
