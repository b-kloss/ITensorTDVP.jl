using ITensors
using ITensorTDVP
using Observers
using HDF5
include("./boson.jl")
include("./model.jl")
include("./gates.jl")


function measure_corr_ij(kets,bras,i,j; conjugate=false)
  if conjugate
    
    return inner(Base.conj(bras[i]),kets[j])  ###lookup inner whether inner automatically applies dag or not
  else
    return inner(bras[i],kets[j])
  end
end



function log(fname, oname, time, data)
  fid=h5open(fname,"a")
  create_dataset(fid,oname+"/t"+string(time),data)
  close(fid)
  return nothing
end

let 
  function myaddition(a,b)
    return a+b
  end  
end




let
tebd_cutoff=1e-18
tebd_dt=0.01
tebd_order = 2  #2 or 4(not implemented yet)
log_dt=0.1
tdvp_dt=0.05
tdvp_cutoff=1e-12
tdvp_cutoff_compress=0.0
tdvp_nsite=1
tdvp_maxdim=64
N = 16
dim=16

phys_bos = siteinds("MyBoson", N,dim=dim,conserve_qns=true,conserve_number=false,)
ancs_bos = siteinds("MyBoson", N,dim=dim, conserve_qns=true,conserve_number=false)
els = siteinds("Fermion",N,conserve_qns=true)
ancs_bos = addtags(ancs_bos,",ancilla")
#sites=Vector{Index{Vector{Pair{QN, Int64}}}}()
sites=Vector{Index}()

for (x,y,z) in zip(phys_bos,els,ancs_bos)
    append!(sites,(x,y,z))
end

global tebd_gates=tfd_holstein_gates(sites, tebd_dt,tebd_order, 1.0,1.0,1.0,0.4)

function tebd_step!(psi::MPS, dummy)
  println("doing TEBD step")
  Nsteps=convert(Int, ceil(abs(log_dt / tebd_dt)))
  for i in range(1,Nsteps)
    for gate in tebd_gates
      psi = apply(gate,psi,cutoff=tebd_cutoff)
    end
  end
  println(maxlinkdim(psi)," ",minimum(linkdims(psi)))
  return psi, dummy
end

function tdvp_step!(ψ::MPS,PH)
  println("doing TDVP step")
  println(typeof(PH))
  ψ = tdvp(
  PH,
  -log_dt*im,
  ψ;
  time_step=-im * tdvp_dt,
  maxdim=tdvp_maxdim,
  cutoff=tdvp_cutoff*(tdvp_nsite==1 ? 1e7 : 1.0),
  cutoff_compress=tdvp_cutoff_compress,
  outputlevel=1,
  nsite=tdvp_nsite,
  (observer!)=obs,
)
  return ψ,PH
end

function propagate!(ψ,PH,propfunc)
  ldims=linkdims(ψ)
  
  #if all([dim>1 for dim in ldims])
  #  ψ=tdvp_step!(ψ,PH)
  ψ,PH=propfunc(ψ,PH)
  #
  #else
  #  ψ=tebd_step!(ψ)
  #  #position!(PH,ψ,1)
  #end
  
  return ψ,PH
end

opsum=tfd_holstein(N;omega=1.0,t=1.0,alpha=1.0, T=0.4, order=["phys_bos","el","anc_bos"])
states = [n == div(N,2) ? [1,"Occ",1] : [1,"Emp",1] for n=1:N]
states=reduce(vcat,states)
@show states
@show typeof(sites[1])
#@show states
#states = [isodd(n) ? "Up" : "Dn" for n=1:N]
#@show typeof(states)
#@show typeof(s)



ψ = MPS(ComplexF64,sites, states)

@show linkdims(ψ)
@show maxlinkdim(ψ)


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
PH=ProjMPO(H)
ψ0=deepcopy(ψ)
for i in range(1,13)
  @time ψ,PH=propagate!(ψ,PH,tebd_step!)
  Nphys_sites=div(length(ψ),3)
  #@show i, expect(ψ, "n"; sites=[(j-1)*3 + 2 for j in 1:Nphys_sites])
  center_value=measure_corr_ij([ψ,],[ψ,],1,1,conjugate=true)
  @show 2*i, center_value

  center_value2=measure_corr_ij([ψ,],[ψ0],1,1,conjugate=false)
  @show i, center_value2
  println(maxlinkdim(ψ)," ",minimum(linkdims(ψ)))
end
#PH=ProjMPO(H)
for i in range(1,10)
  if maxlinkdim(ψ)==tdvp_maxdim
    tdvp_cutoff_compress=1e-12
  end
  @time ψ,H=propagate!(ψ,H,tdvp_step!)
end

@show inner(ϕ', H, ϕ) / inner(ϕ, ϕ)
res = results(obs)
steps = res["steps"]
times = res["times"]
psis = res["pops"]
maxdim = res["maxdim"]

@show res
end