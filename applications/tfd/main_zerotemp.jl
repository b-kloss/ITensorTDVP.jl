using ITensors
using ITensorTDVP
using Observers
include("./boson.jl")
N = 10
dim=8
s = siteinds("S=1/2", N, conserve_qns=true)
#@show s,typeof(s)
phys_bos = siteinds("MyBoson", N,dim=8,conserve_qns=false,conserve_number=false,)
els = siteinds("Fermion",N,conserve_qns=false)
#sites=Vector{Index{Vector{Pair{QN, Int64}}}}()
sites=Vector{Index}()

for (x,y) in zip(phys_bos,els)
    append!(sites,(x,y))
end

states = [n == div(N,2) ? [1,"Occ"] : [1,"Emp"] for n=1:N]

states=reduce(vcat,states)
@show states

    

function holstein(n;omega=1.0,t=1.0, alpha=1.0, T=0.0, order=["phys_bos","el","anc_bos"])
    os = OpSum()
    @assert order==["phys_bos","el"]
    elpos=2
    ancpos=3
    physpos=1
    for j in 1:(n - 1)
      elj=elpos
      os += t,"C", 2*(j-1)+elj, "Cdag", 2*(j) +elj
      os += t,"Cdag", 2*(j-1)+elj, "C", 2*(j) +elj
    end
    for j in 1:N
      os += alpha,"A", 2*(j-1) + physpos, "N", 2*(j-1)+elpos   #V
      os += alpha, "Adag", 2*(j-1) + physpos, "N", 2*(j-1)+elpos
      os += omega, "N", 2*(j-1)+physpos #local oscillator
    end
    return os
  end

ψ = MPS(ComplexF64,sites, states)
ψ2 = randomMPS(ComplexF64,sites, states;linkdims=10)

ψ = ψ + 1e-4 * ψ2
Nphys_sites=div(length(ψ),2)
@show  expect(ψ, "N"; sites=[(2*(j-1) + 2) for j in 1:Nphys_sites])
opsum=holstein(n,omega=1.0,t=1.0,alpha=0.0, T=0.4, order=["phys_bos","el"]),sites
#@show typeof(opsum[2])
H = MPO(opsum[1],sites)
H2 = MPO(opsum[1],opsum[2])
#@show siteinds(H)
#@show siteinds(H2)

function measure_pop(; psi, bond, half_sweep)
    if bond == 1 && half_sweep == 2
       Nphys_sites=div(length(psi),2)
      return expect(psi, "n"; sites=[(j-1)*2 + 2 for j in 1:Nphys_sites])
    end
    return nothing
end

function bonddim(;psi,bond,half_sweep)
    if bond == 1 && half_sweep == 2
       return maxlinkdim(psi)
    end
    return nothing
end
function step(; sweep, bond, half_sweep)
    if bond == 1 && half_sweep == 2
      return sweep
    end
    return nothing
  end
function current_time(; current_time, bond, half_sweep)
    if bond == 1 && half_sweep == 2
      return current_time
    end
    return nothing
  end
obs = Observer(
    "pops" => measure_pop,  "steps" => step, "times" => current_time,"maxdim" => bonddim
  )
@show inner(ψ', H2, ψ) / inner(ψ, ψ)
cutoff = 1e-12
tau = 0.1
ttotal = 2.0

ϕ = tdvp(
  H,
  -ttotal*im,
  ψ;
  time_step=-im * tau,
  maxdim=128,
  cutoff=cutoff,
  outputlevel=1,
  nsite=2,
  (observer!)=obs,
)

@show inner(ϕ', H, ϕ) / inner(ϕ, ϕ)
res = results(obs)
steps = res["steps"]
times = res["times"]
psis = res["pops"]
maxdim = res["maxdim"]

@show res
