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

# FLEX approximation
Author: [Kosuke Nogaki](mailto:nogaki.kosuke.83v@st.kyoto-u.ac.jp)

## Theory of FLEX in the paramagnetic state

+++

## Code implementation

```{code-cell}
using PyPlot
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["text.usetex"] = true
rcParams["font.family"] = "serif"
rcParams["font.size"] = 16
rcParams["text.latex.preamble"] = raw"\usepackage{amsmath}"

using Revise
using FFTW
using LinearAlgebra
using Roots
using SparseIR
import SparseIR: Statistics
#using JET

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
```

### Parameter setting

```{code-cell}
### System parameters
t    = 1      # hopping amplitude
W    = 8*t    # bandwidth
wmax = 10     # set wmax >= W

T    = 0.1    # temperature
beta = 1/T    # inverse temperature
n    = 0.85   # electron filling, here per spin per lattice site (n=1: half filling)
U    = 4.0    # Hubbard interaction

### Numerical parameters
nk1, nk2  = 24, 24    # number of k_points along one repiprocal crystal lattice direction k1 = kx, k2 = ky
nk        = nk1*nk2
IR_tol    = 1e-10     # accuary for l-cutoff of IR basis functions
sfc_tol   = 1e-4      # accuracy for self-consistent iteration
maxiter   = 30        # maximal number of iterations in self-consistent cycle
mix       = 0.2       # mixing parameter for new 
U_maxiter = 50        # maximal number of iteration steps in U renormalization loop
;
```

### Generating meshes

```{code-cell}
"""
Holding struct for k-mesh and sparsely sampled imaginary time 'tau' / Matsubara frequency 'iw_n' grids.
Additionally we defines the Fourier transform routines 'r <-> k'  and 'tau <-> l <-> wn'.
 """
struct Mesh
    nk1         ::Int64
    nk2         ::Int64
    nk          ::Int64
    ek          ::Array{Float64,2} 
    iw0_f       ::Int64
    iw0_b       ::Int64
    fnw         ::Int64
    fntau       ::Int64
    bnw         ::Int64
    bntau       ::Int64
    IR_basis_set::FiniteTempBasisSet
end

"""Initiarize function"""
function Mesh(
        nk1         ::Int64,
        nk2         ::Int64,
        IR_basis_set::FiniteTempBasisSet,
        )::Mesh
    
    nk::Int64 = nk1*nk2

    # Compute Hamiltonian
    ek = Array{ComplexF64,2}(undef, nk1, nk2)
    for iy in 1:nk2, ix in 1:nk1
        kx::Float64 = (2*π*(ix-1))/nk1
        ky::Float64 = (2*π*(iy-1))/nk2
        ek[ix, iy] = -2.0*(cos(kx)+cos(ky)) 
    end
    
    # lowest Matsubara frequency index
    iw0_f = findall(x->x==1,IR_basis_set.smpl_wn_f.sampling_points)[1]
    iw0_b = findall(x->x==0,IR_basis_set.smpl_wn_b.sampling_points)[1]
    
    # the number of sampling point for fermion and boson
    fnw   = length(IR_basis_set.smpl_wn_f.sampling_points)
    fntau = length(IR_basis_set.smpl_tau_f.sampling_points)
    bnw   = length(IR_basis_set.smpl_wn_b.sampling_points)
    bntau = length(IR_basis_set.smpl_tau_b.sampling_points)

    # Return
    Mesh(nk1, nk2, nk, ek, iw0_f, iw0_b, fnw, fntau, bnw, bntau, IR_basis_set)
end

function smpl_obj(mesh::Mesh, statistics::SparseIR.Statistics)
    """ Return sampling object for given statistic """
    if statistics == fermion
        smpl_tau = mesh.IR_basis_set.smpl_tau_f
        smpl_wn  = mesh.IR_basis_set.smpl_wn_f
    elseif statistics == boson
        smpl_tau = mesh.IR_basis_set.smpl_tau_b
        smpl_wn  = mesh.IR_basis_set.smpl_wn_b
    end    
    return smpl_tau, smpl_wn
end

"""Fourier transformation"""    
function tau_to_wn(mesh::Mesh, statistics, obj_tau) where {T <: SparseIR.Statistics}
    """ Fourier transform from tau to iw_n via IR basis """
    smpl_tau, smpl_wn = smpl_obj(mesh, statistics)

    obj_l = fit(smpl_tau, obj_tau, dim=1)
    obj_wn = evaluate(smpl_wn, obj_l, dim=1)
    return obj_wn
end

function wn_to_tau(mesh::Mesh, statistics::Statistics, obj_wn)
    """ Fourier transform from iw_n to tau via IR basis """
    smpl_tau, smpl_wn = smpl_obj(mesh, statistics)

    obj_l   = fit(smpl_wn, obj_wn, dim=1)
    obj_tau = evaluate(smpl_tau, obj_l, dim=1)
    return obj_tau
end
 
function k_to_r(mesh::Mesh, obj_k)
    """ Fourier transform from k-space to real space """
    obj_r = fft(obj_k,[2,3])
    return obj_r
end

function r_to_k(mesh::Mesh, obj_r)
    """ Fourier transform from real space to k-space """
    obj_k = ifft(obj_r,[2,3])/mesh.nk
    return obj_k
end

@assert typestable(tau_to_wn, (Mesh, SparseIR.Statistics, Array{ComplexF64,4}))
@assert typestable(wn_to_tau, (Mesh, SparseIR.Statistics, Array{ComplexF64,4}))
```

### FLEX loop solver

```{code-cell}
"""
Solver struct to calculate the FLEX loop self-consistently.
After initializing the Solver by `solver = FLEXSolver(mesh, beta, U, n, sigma_init, sfc_tol, maxiter, U_maxiter, mix)'
it can be run by `solve(solver)`.
 """
mutable struct FLEXSolver
    mesh     ::Mesh
    beta     ::Float64
    U        ::Float64
    n        ::Float64
    sfc_tol  ::Float64
    maxiter  ::Int64
    U_maxiter::Int64
    mix      ::Float64
    verbose  ::Bool
    mu       ::Float64
    gkio     ::Array{ComplexF64,3}
    grit     ::Array{ComplexF64,3}
    ckio     ::Array{ComplexF64,3}
    V        ::Array{ComplexF64,3}
    sigma    ::Array{ComplexF64,3}
end

"""Initiarize function"""
function FLEXSolver(
        mesh      ::Mesh,
        beta      ::Float64,
        U         ::Float64,
        n         ::Float64,
        sigma_init::Array{ComplexF64,3};
        sfc_tol   ::Float64=1e-4,
        maxiter   ::Int64  =100,
        U_maxiter ::Int64  =10,
        mix       ::Float64=0.2,
        verbose   ::Bool   =true
        )::FLEXSolver
    
        mu::Float64 = 0.0
    
        gkio  = Array{ComplexF64}(undef, mesh.fnw,   mesh.nk1, mesh.nk2)
        grit  = Array{ComplexF64}(undef, mesh.fntau, mesh.nk1, mesh.nk2)
        ckio  = Array{ComplexF64}(undef, mesh.bnw,   mesh.nk1, mesh.nk2)
        V     = Array{ComplexF64}(undef, mesh.bntau, mesh.nk1, mesh.nk2)
        sigma = sigma_init
    
        solver = FLEXSolver(mesh, beta, U, n, sfc_tol, maxiter, U_maxiter, mix, verbose, mu, gkio, grit, ckio, V, sigma)
    
        solver.mu = mu_calc(solver)
        gkio_calc(solver,solver.mu)
        grit_calc(solver)
        ckio_calc(solver)
        return solver
end

#%%%%%%%%%%% Loop solving instance
function solve(solver::FLEXSolver)
    """ FLEXSolver.solve() executes FLEX loop until convergence """
    # check whether U < U_crit! Otherwise, U needs to be renormalized.
    if maximum(abs, solver.ckio) * solver.U >= 1
        U_renormalization(solver)
    end
            
    # perform loop until convergence is reached:
    for it in 1:solver.maxiter
        sigma_old = copy(solver.sigma)
        loop(solver)
        
        # check whether solution is converged.
        sfc_check = sum(abs.(solver.sigma-sigma_old))/sum(abs.(solver.sigma))

        if solver.verbose
            println(it, '\t', sfc_check)
        end
        if sfc_check < solver.sfc_tol
            println("FLEX loop converged at desired accuracy")
            break
        end
    end
end
    
function loop(solver::FLEXSolver)
    """ FLEX loop """
    #t1 = time_ns()
    gkio_old = copy(solver.gkio)
    
    #t2 = time_ns()

    V_calc(solver)
    sigma_calc(solver)
    
    #t3 = time_ns()
        
    solver.mu = mu_calc(solver)
    #t4 = time_ns()
    gkio_calc(solver,solver.mu)
    #t5 = time_ns()
    
    solver.gkio .= solver.mix*solver.gkio .+ (1-solver.mix)*gkio_old
        
    grit_calc(solver)
    ckio_calc(solver)
    #t6 = time_ns()
    #println("debug ", (t2-t1)*1e-9, " ", (t3-t2)*1e-9, (t4-t3)*1e-9, " ", (t5-t4)*1e-9, " ", (t6-t5)*1e-9)
end


#%%%%%%%%%%% U renormalization loop instance
function U_renormalization(solver::FLEXSolver)
    """ Loop for renormalizing U if Stoner enhancement U*max{chi0} >= 1. """
    println("WARNING: U is too large and the spin susceptibility denominator will diverge/turn unphysical!")
    println("Initiate U renormalization loop.")
    
    # save old U for later
    U_old::Float64 = solver.U
    # renormalization loop may run infinitely! Insert break condition after U_it_max steps
    U_it::Int64 = 0
    
    while U_old*maximum(abs, solver.ckio) >= 1.0
        U_it += 1
        
        # remormalize U such that U*chi0 < 1
        solver.U = solver.U / (maximum(abs, solver.ckio)*solver.U + 0.01)
        println(U_it, '\t', solver.U, '\t', U_old)
        
        # perform one shot FLEX loop
        loop(solver)
        
        # reset U
        solver.U = U_old
        
        # break condition for too many steps
        if U_it == solver.U_maxiter
            println("U renormalization reached breaking point")
            break
        end
    end
    println("Leaving U renormalization...")
end

#%%%%%%%%%%% Calculation steps
function gkio_calc(solver::FLEXSolver, mu::Float64)
    """ calculate Green function G(iw,k) """
    for iy in 1:solver.mesh.nk2, ix in 1:solver.mesh.nk1, iw in 1:solver.mesh.fnw
        iv::ComplexF64 = (im * π/solver.beta) * solver.mesh.IR_basis_set.smpl_wn_f.sampling_points[iw]
        solver.gkio[iw,ix,iy] = 1.0/(iv - solver.mesh.ek[ix, iy] + mu - solver.sigma[iw,ix,iy])
    end
end

function grit_calc(solver::FLEXSolver)
    """ Calculate real space Green function G(tau,r) [for calculating chi0 and sigma] """
    # Fourier transform
    grio = k_to_r(solver.mesh, solver.gkio)
    solver.grit .= wn_to_tau(solver.mesh, fermion, grio)
end

function ckio_calc(solver::FLEXSolver)
    """ Calculate irreducible susciptibility chi0(iv,q) """
    crit = Array{ComplexF64}(undef, solver.mesh.bntau, solver.mesh.nk1, solver.mesh.nk2)
    for iy in 1:solver.mesh.nk2, ix in 1:solver.mesh.nk1, it in 1:solver.mesh.bntau
        crit[it,ix,iy] = solver.grit[it,ix,iy] * solver.grit[solver.mesh.bntau-it+1,ix,iy]
    end

    # Fourier transform
    ckit = r_to_k(solver.mesh, crit)
    solver.ckio .= tau_to_wn(solver.mesh, boson, ckit)
end

function V_calc(solver::FLEXSolver)
    """ Calculate interaction V(tau,r) from RPA-like spin and charge susceptibility for calculating sigma """
    # check whether U is too large and give warning
    if maximum(abs.(solver.ckio))*solver.U >= 1
        error("U*max(chi0) >= 1! Paramagnetic phase is left and calculations will turn unstable!")
    end
        
    # spin and charge susceptibility
    chi_spin   = solver.ckio ./ (1 .- solver.U .* solver.ckio)
    chi_charge = solver.ckio ./ (1 .+ solver.U .* solver.ckio)

    Vkio = (1.5*solver.U^2) .* chi_spin .+ (0.5*solver.U^2) .* chi_charge .- (solver.U^2) .* solver.ckio
    # Constant Hartree Term V ~ U needs to be treated extra, since they cannot be modeled by the IR basis.
    # In the single-band case, the Hartree term can be absorbed into the chemical potential.

    # Fourier transform
    Vrio = k_to_r(solver.mesh, Vkio)
    solver.V .= wn_to_tau(solver.mesh, boson, Vrio)
end

function sigma_calc(solver::FLEXSolver)
    """ Calculate self-energy Sigma(iw,k) """
    sigmarit = solver.V .* solver.grit
    
    # Fourier transform
    sigmakit = r_to_k(solver.mesh, sigmarit)
    solver.sigma .= tau_to_wn(solver.mesh, fermion, sigmakit)
end
    
    
#%%%%%%%%%%% Setting chemical potential mu
function calc_electron_density(solver::FLEXSolver,mu::Float64)::Float64
    """ Calculate electron density from Green function """
    #t1 = time_ns()
    gkio_calc(solver,mu)
    gio = dropdims(sum(solver.gkio,dims=(2,3)),dims=(2,3))/solver.mesh.nk
    #t2 = time_ns()
    
    g_l = fit(solver.mesh.IR_basis_set.smpl_wn_f,gio, dim=1)
    #t3 = time_ns()
    g_tau0 = dot(solver.mesh.IR_basis_set.basis_f.u(0), g_l)
    #t4 = time_ns()
    #println("timing ", 1e-9*(t2-t1),  " ", 1e-9*(t3-t2), " ", 1e-9*(t4-t3))
    
    n  = 1.0 + real(g_tau0)
    n  = 2.0 * n #for spin
end

function mu_calc(solver::FLEXSolver)::Float64
    """ Find chemical potential for a given filling n0 via brent's root finding algorithm """
    f  = x -> calc_electron_density(solver,x) - solver.n
        
    mu = find_zero(f, (3*minimum(solver.mesh.ek), 3*maximum(solver.mesh.ek)), Roots.Brent()) 
end
@assert typestable(U_renormalization, FLEXSolver)
@assert typestable(solve, FLEXSolver)
```

### Execute FLEX loop

```{code-cell}
# initialize calculation
#t1 = time_ns()
IR_basis_set = FiniteTempBasisSet(beta, wmax, IR_tol)
#t2 = time_ns()
mesh = Mesh(nk1, nk2, IR_basis_set)
sigma_init = zeros(ComplexF64,(mesh.fnw, nk1, nk2))
solver = FLEXSolver(mesh, beta, U, n, sigma_init, sfc_tol=sfc_tol, maxiter=maxiter, U_maxiter=U_maxiter, mix=mix)
#t3 = time_ns()

# perform FLEX loop
solve(solver)
#t4 = time_ns()
#println((t2-t1)*1e-9)
#println((t3-t2)*1e-9)
#println((t4-t3)*1e-9)
```

```{code-cell}
#using BenchmarkTools
#@benchmark loop(solver)
```

```{code-cell}
#@benchmark gkio_calc(solver, 0.0)
#@code_warntype gkio_calc(solver, 0.0)
```

```{code-cell}
#@benchmark calc_electron_density(solver, 0.0)
```

```{code-cell}
#@report_opt target_modules=(@__MODULE__,)  U_renormalization(solver)
#using BenchmarkTools
#@benchmark loop(solver)
#@benchmark V_calc(solver)
```

```{code-cell}
#@benchmark sigma_calc(solver)
```

```{code-cell}
#@benchmark mu_calc(solver)
```

#### Visualize results

```{code-cell}
# plot 2D k-dependence of lowest Matsubara frequency of e.g. green function
myx = (2 .* collect(1:nk1) .- 1) ./ nk1
myy = (2 .* collect(1:nk1) .- 1) ./ nk2
plt.figure(figsize=(3.7,3))
plt.pcolormesh(myx, myy, real.(solver.gkio[solver.mesh.iw0_f,:,:]'), shading="auto")
ax = plt.gca()
ax.set_xlabel(L"k_x/\pi")
ax.set_xlim([0,2])
ax.set_ylabel(L"k_y/\pi")
ax.set_ylim([0,2])
ax.set_aspect("equal")
plt.colorbar()
plt.show()
```

```{code-cell}
# plot 2D k-dependence of lowest Matsubara frequency of e.g. chi0
plt.figure(figsize=(3.7,3))
plt.pcolormesh(myx, myy, real.(solver.ckio[solver.mesh.iw0_b,:,:]'), shading="auto")
ax = plt.gca()
ax.set_xlabel(L"k_x/\pi")
ax.set_xlim([0,2])
ax.set_ylabel(L"k_y/\pi")
ax.set_ylim([0,2])
ax.set_aspect("equal")
plt.colorbar()
plt.show()
```

## Linearized Eliashberg equation

+++

### Code implementation

#### Linearized Eliashberg solver

```{code-cell}
"""
Solver struct for solving the linearized gap equation using the power method.
It takes FLEX results as an input.
"""
mutable struct LinearizedGapSolver
    mesh      ::Mesh
    gkio      ::Array{ComplexF64,3}
    V_singlet ::Array{ComplexF64,3}
    delta     ::Array{ComplexF64,3}
    frit      ::Array{ComplexF64,3}
    U         ::Float64
    maxiter   ::Int64
    sfc_tol   ::Float64
    verbose   ::Bool
    lam       ::Float64
end

function LinearizedGapSolver(
        FLEX_solver::FLEXSolver;
        maxiter    ::Int64  =50,
        sfc_tol    ::Float64=1e-4,
        verbose    ::Bool   =true
        )::LinearizedGapSolver
        
    ## Initialize necessary quantities from FLEX loop
    mesh       = FLEX_solver.mesh
    gkio       = FLEX_solver.gkio
    U          = FLEX_solver.U
        
    maxiter = maxiter
    sfc_tol = sfc_tol
    verbose = verbose
        
    ## Initialize trial gap function
    # Here we focus on a d-wave symmetric solution
    delta = Array{ComplexF64}(undef, FLEX_solver.mesh.fnw,   FLEX_solver.mesh.nk1, FLEX_solver.mesh.nk2)
    frit  = Array{ComplexF64}(undef, FLEX_solver.mesh.fntau, FLEX_solver.mesh.nk1, FLEX_solver.mesh.nk2)
    for iy in 1:FLEX_solver.mesh.nk2, ix in 1:FLEX_solver.mesh.nk1, iw in 1:FLEX_solver.mesh.fnw
        kx::Float64 = (2*π*(ix-1))/FLEX_solver.mesh.nk1
        ky::Float64 = (2*π*(iy-1))/FLEX_solver.mesh.nk2
        delta[iw,ix,iy] = cos(kx) - cos(ky) 
    end
    
    #normalize initial guess
    normalize!(delta)   
        
    # Initialize interaction
    V_singlet = V_singlet_calc(FLEX_solver)
        
    ## Initialize eigenvalue
    lam::Float64 = 0.0
    gap_solver = LinearizedGapSolver(mesh, gkio, V_singlet, delta, frit, U, maxiter, sfc_tol, verbose, lam)
end

function solve(gap_solver::LinearizedGapSolver)
    """ Solving instance to find eigenvalue from power method """
    for it in 1:gap_solver.maxiter
        lam_old = gap_solver.lam
        delta_old = copy(gap_solver.delta)
    
        frit_calc(gap_solver)
        deltarit = gap_solver.V_singlet .* gap_solver.frit
    
        # Fourier transform to momentum space
        deltakit = r_to_k(gap_solver.mesh, deltarit)
        gap_solver.delta .= tau_to_wn(gap_solver.mesh, fermion, deltakit)
    
        # calculate eigenvalue
        gap_solver.lam = sum(real.(conj.(gap_solver.delta).* delta_old))
        
        normalize!(gap_solver.delta) 
    
        if gap_solver.verbose
            println(it, '\t', gap_solver.lam, '\t', abs(gap_solver.lam-lam_old))
        end
        if abs(gap_solver.lam-lam_old) < gap_solver.sfc_tol
            break   
        end
    end
end

    
#%%%%%%%%%%% Calculation steps
function V_singlet_calc(solver::FLEXSolver)::Array{ComplexF64,3}
    """ Set up interaction in real space and imaginary time """
    chi_spin   = solver.ckio ./ (1 .- solver.U*solver.ckio)
    chi_charge = solver.ckio ./ (1 .+ solver.U*solver.ckio)
    
    Vkio = 1.5*solver.U^2 * chi_spin .- 0.5*solver.U^2 * chi_charge
    # Constant Hartree Term V ~ U needs to be treated extra, since they cannot be modeled by the IR basis.
    # In the special case of d-wave symmetry, it can be neglected.
    
    # Fourier transform
    Vrio = k_to_r(solver.mesh, Vkio)
    V_singlet = wn_to_tau(solver.mesh, boson, Vrio)
    return V_singlet
end

function frit_calc(gap_solver::LinearizedGapSolver)
    """ Calculate (linearized) anomalous Green function F = |G|^2 * delta for evaluating the gap equation """
    fkio = - gap_solver.gkio.*conj(gap_solver.gkio).*gap_solver.delta
        
    # Fourier transform
    frit = k_to_r(gap_solver.mesh, fkio)
    gap_solver.frit = wn_to_tau(gap_solver.mesh, fermion, frit)
end
```

#### Executing the gap equation solver

```{code-cell}
gap_solver = LinearizedGapSolver(solver, maxiter=maxiter, sfc_tol=sfc_tol)
solve(gap_solver)
println("The superconducting eigenvalue at T=",T," is lambda_d=",gap_solver.lam)
```

#### Visualize results

```{code-cell}
# plot 2D k-dependence of lowest Matsubara frequency of the gap vs. initial guess
delta0 = Array{ComplexF64}(undef, gap_solver.mesh.nk1, gap_solver.mesh.nk2)
for iy in 1:gap_solver.mesh.nk2, ix in 1:gap_solver.mesh.nk1
    kx::Float64 = (2*π*(ix-1))/nk1
    ky::Float64 = (2*π*(iy-1))/nk2
    delta0[ix,iy] = cos(kx) - cos(ky) 
end
normalize!(delta0)

plt.figure(figsize=(3.7,3))
plt.pcolormesh(myx, myy, real.(delta0[:,:]'), cmap="RdBu", shading="auto")
ax = plt.gca()
#plt.pcolormesh(myx, myy, real.(gap_solver.delta[gap_solver.mesh.iw0_f,:,:]'))
ax.set_xlabel(L"k_x/\pi")
ax.set_xlim([0,2])
ax.set_ylabel(L"k_y/\pi")
ax.set_ylim([0,2])
ax.set_aspect("equal")
ax.set_title(L"\Delta^0_d(k)")
plt.colorbar()


plt.figure(figsize=(3.7,3))
plt.pcolormesh(myx, myy, real.(gap_solver.delta[gap_solver.mesh.iw0_f,:,:]'), cmap="RdBu", shading="auto")
ax = plt.gca()
ax.set_xlabel(L"k_x/\pi")
ax.set_xlim([0,2])
ax.set_ylabel(L"k_y/\pi")
ax.set_ylim([0,2])
ax.set_aspect("equal")
ax.set_title(L"\Delta_d(k)")
plt.colorbar()
plt.show()
```

```{code-cell}
# plot 2D k-dependence of lowest Matsubara frequency of the anomalous Green function
fkio = - gap_solver.gkio.*conj(gap_solver.gkio).*gap_solver.delta
plt.figure(figsize=(3.7,3))
plt.pcolormesh(myx, myy, real.(fkio[mesh.iw0_f,:,:]'), cmap="RdBu", shading="auto")
ax = plt.gca()
ax.set_xlabel(L"k_x/\pi")
ax.set_xlim([0,2])
ax.set_ylabel(L"k_y/\pi")
ax.set_ylim([0,2])
ax.set_aspect("equal")
ax.set_title(L"\Delta_d(k)")
plt.colorbar()
plt.show()
```

## Example: Antiferromagnetic fluctuations and $d$-wave superconductivity in the square-lattice Hubbard model

```{code-cell}
#%%%%%%%%%%%%%%% Parameter settings
println("Initialization...")
# system parameters
t    = 1      # hopping amplitude
n    = 0.85   # electron filling, here per spin per lattice site (n=1: half filling)
U    = 4.0    # Hubbard interaction

W    = 8*t    # bandwidth
wmax = 10     # set wmax >= W
T_values = [0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.025]# temperature

# numerical parameters
nk1, nk2  = 64, 64    # k-mesh sufficiently dense!
nk        = nk1*nk2
IR_Lambda = 10^3      # dimensionless IR parameter >= w_max * beta_min = 400
IR_tol    = 1e-7      # accuary for l-cutoff of IR basis functions
sfc_tol   = 1e-4      # accuracy for self-consistent iteration
it_max    = 30        # maximal number of iterations in self-consistent cycle
mix       = 0.2       # mixing parameter for new 
U_it_max  = 50        # maximal number of iteration steps in U renormalization loop

# initialize first IR basis set (no recalculation afterwrds)
beta_init = 1.0/T_values[1]
IR_basis_set = FiniteTempBasisSet(beta_init, IR_Lambda/beta_init, IR_tol)

# set initial self_energy - will be set to previous calculation step afterwards
sigma_init = zeros(ComplexF64,(length(IR_basis_set.smpl_wn_f.sampling_points), nk1, nk2))

# empty arrays for results
lam_T     = Array{Float64}(undef,length(T_values))
chiSmax_T = Array{Float64}(undef,length(T_values))
chi_s_plt = Array{Float64}(undef,nk1,nk2)

#%%%%%%%%%%%%%%% Calculations for different T values
for T_it in 1:length(T_values)
    T = T_values[T_it]
    println("Now: T = ", T)
    beta = 1/T

    # initialize meshes
    IR_basis_set = FiniteTempBasisSet(beta, IR_Lambda/beta, IR_tol, sve_result=IR_basis_set.basis_f.sve_result)
    mesh = Mesh(nk1, nk2, IR_basis_set)
    
    # calculate FLEX loop
    solver = FLEXSolver(mesh, beta, U, n, sigma_init, sfc_tol=sfc_tol, maxiter=it_max, U_maxiter=U_it_max, mix=mix, verbose=false)
    solve(solver)
    
    sigma_init = copy(solver.sigma)
    
    # calculate linearized gap equation
    gap_solver = LinearizedGapSolver(solver, maxiter=it_max, sfc_tol=sfc_tol, verbose=false)
    solve(gap_solver)
    
    # save data for plotting
    lam_T[T_it] = gap_solver.lam
    chi_spin   = solver.ckio ./ (1 .- solver.U*solver.ckio)
    chiSmax_T[T_it] = real(findmax(abs.(chi_spin))[1])
    
    if T == 0.03
        chi_s_plt .= real.(chi_spin[solver.mesh.iw0_b,:,:])
    end

end
```

```{code-cell}
fig   = plt.figure(figsize=(10,4),constrained_layout=true)
spec  = PyPlot.matplotlib.gridspec.GridSpec(ncols=2, nrows=1, figure=fig)
f_ax1 = fig.add_subplot(spec[0, 1])
f_ax2 = fig.add_subplot(spec[0, 0])

# first panel with momentum dependence of static spin susceptibility
chi_s_HSP= Array{Float64}(undef,Int64(3*(solver.mesh.nk1/2)))
for i in 1:Int64(solver.mesh.nk1/2)
    chi_s_HSP[i] = chi_s_plt[1,i]
    chi_s_HSP[Int64(solver.mesh.nk1/2)+i] = chi_s_plt[i,Int64(solver.mesh.nk1/2)+1]
    chi_s_HSP[solver.mesh.nk1+i] = chi_s_plt[Int64(solver.mesh.nk1/2)-i+2,Int64(solver.mesh.nk1/2)-i+2]
end

k_HSP= Array{Float64}(undef,Int64(3*(solver.mesh.nk1/2)))
for i in 1:solver.mesh.nk1
    k_HSP[i] = 2*(i-1)/nk1
end
for i in 1:Int64(solver.mesh.nk1/2)
    k_HSP[solver.mesh.nk1+i] = 2*(solver.mesh.nk1-1)/solver.mesh.nk1 + sqrt(2)*2*(i-1)/solver.mesh.nk1
end

f_ax1.plot(k_HSP, chi_s_HSP,"-")
f_ax1.set_xlim([0,2+sqrt(2)])
f_ax1.set_xticks([0,1,2,2+sqrt(2)])
f_ax1.set_xticklabels([L"\Gamma",L"X",L"M",L"\Gamma"])
f_ax1.set_ylim([0,26])
f_ax1.set_xlabel("")
f_ax1.set_ylabel(L"\chi_{\mathrm{s}}(i\nu=0,{\bf{q}})", fontsize=14)
f_ax1.grid()

Ones      = ones(Int,length(T_values))
# second panel with T-dependence of lambda_d and 1/chi_s,max
f_ax2.plot(T_values, lam_T, "-x", label=L"\lambda_d")
f_ax2.plot(T_values, Ones./chiSmax_T, marker="x", label=L"1/\chi_{\mathrm{s},\mathrm{max}}")
f_ax2.set_xlim([0.01,0.08])
f_ax2.set_ylim([0,1])
f_ax2.set_xlabel(L"T/t", fontsize=14)
f_ax2.set_ylabel(L"\lambda_d, 1/\chi_{\mathrm{s},\mathrm{max}}", fontsize=14)
f_ax2.legend()
f_ax2.grid()
plt.show()
```

```{code-cell}

```

```{code-cell}

```
