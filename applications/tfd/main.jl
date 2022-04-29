using ITensors
using ITensorTDVP
using Observers
include("./boson.jl")
include("./model.jl")
include("./gates.jl")

let

N = 16
dim=16
s = siteinds("S=1/2", N, conserve_qns=true)
#@show s,typeof(s)
phys_bos = siteinds("MyBoson", N,dim=dim,conserve_qns=true,conserve_number=false,)
ancs_bos = siteinds("MyBoson", N,dim=dim, conserve_qns=true,conserve_number=false)
els = siteinds("Fermion",N,conserve_qns=true)
ancs_bos = addtags(ancs_bos,",ancilla")
#sites=Vector{Index{Vector{Pair{QN, Int64}}}}()
sites=Vector{Index}()

for (x,y,z) in zip(phys_bos,els,ancs_bos)
    append!(sites,(x,y,z))
end

gates=tfd_holstein_gates(sites, 0.01,2, 1.0,1.0,1.0,0.4)
states = [n == div(N,2) ? [1,"Occ",1] : [1,"Emp",1] for n=1:N]
N=N
opsum=tfd_holstein(N;omega=1.0,t=1.0,alpha=1.0, T=0.4, order=["phys_bos","el","anc_bos"])
states=reduce(vcat,states)
@show states
@show typeof(sites[1])
#@show states
#states = [isodd(n) ? "Up" : "Dn" for n=1:N]
#@show typeof(states)
#@show typeof(s)



ψ = MPS(ComplexF64,sites, states)
for i in 1:100
  for gate in gates
    ψ = apply(gate,ψ,cutoff=1e-16)
  end
end
@show linkdims(ψ) 
#ψ2 = randomMPS(ComplexF64,sites, states;linkdims=16)
#normalize!(ψ2)
#ψ = ψ + 1e-8 * ψ2
#ψ2 = randomMPS(ComplexF64,sites;linkdims=8)
#normalize!(ψ2)
#ψ = ψ + 1e-8 * ψ2
@show maxlinkdim(ψ)
#@show ψ[4]
#@show ψ[5]
#@show ψ[6]


Nphys_sites=div(length(ψ),3)
@show  expect(ψ, "N"; sites=[(3*(j-1) + 2) for j in 1:Nphys_sites])
@show typeof(opsum)
#@show typeof(opsum[2])
H = MPO(opsum,sites)
#@show siteinds(H)
#@show siteinds(H2)

function measure_pop(; psi, bond, half_sweep)
    if bond == 1 && half_sweep == 2
       Nphys_sites=div(length(psi),3)
      return expect(psi, "n"; sites=[(j-1)*3 + 2 for j in 1:Nphys_sites])
    end
    return nothing
end

function bonddim(;psi,bond,half_sweep)
    if bond == 1 && half_sweep == 2
       return linkdims(psi)
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
@show inner(ψ', H, ψ) / inner(ψ, ψ)
cutoff = 1e-12
tau = 0.1
ttotal = 2.0

ϕ = tdvp(
  H,
  -ttotal*im,
  ψ;
  time_step=-im * tau,
  maxdim=64,
  cutoff=cutoff*10e6,
  cutoff_compress=cutoff,
  outputlevel=1,
  nsite=1,
  (observer!)=obs,
)

@show inner(ϕ', H, ϕ) / inner(ϕ, ϕ)
res = results(obs)
steps = res["steps"]
times = res["times"]
psis = res["pops"]
maxdim = res["maxdim"]

@show res
end