using ITensors
using ITensorTDVP
using HDF5
import MPI
MPI.Init()
comm = MPI.COMM_WORLD


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

function measure_corr(ket,bra;conjugate=true)
    if conjugate
        return inner(Base.conj(bra),ket)  ###lookup inner whether inner automatically applies dag or not
    else
        return inner(bra,ket)
    end
 end

 function measure_corr(i::Int,sendpsi::MPS,recvpsi::MPS,sendpsi0::MPS)
    root=i
    recvpsi=MPI.bcast(sendpsi,root,comm)
    gf = measure_corr(sendpsi,recvpsi,conjugate=true)
    gfhalf = measure_corr(sendpsi0,recvpsi,conjugate=false)
    allgf=MPI.Gather(gf, root,comm)
    allgfhalf=MPI.Gather(gfhalf, root,comm)
    return allgf,allgfhalf
end

function measure_corr(sendpsi::MPS,recvpsi::MPS,sendpsi0::MPS)
    Nphys_sites=div(length(sendpsi),3)
    results=Vector{Vector{ComplexF64}()}
    results_half=Vector{Vector{ComplexF64}()}
    for i in 0:Nphys_sites-1
        gfs,gfs_half=measure_corr(i,sendpsi,recvpsi,sendpsi0)
        push!(results,gfs)
        push!(results_half,gfs_half)
    end
    return results,results_half
end


function log(fname, oname, time, data)
  fid=h5open(fname,"r+")
  write(fid,oname*"/t"*string(time),data)
  close(fid)
  return nothing
end





let
outf="data.h5"
outfid=h5open("data.h5","w")
close(outfid)
tebd_cutoff=1e-14
tebd_dt=0.05
tebd_order = 2  #2 or 4(4 doesn't seem to work as well as it should)
log_dt=0.1
tdvp_dt=0.05
tdvp_cutoff=1e-12
tdvp_cutoff_compress=0.0
tdvp_nsite=2
tdvp_maxdim=64
final_time = 1.0
N = 4
dim=16

total_Nsteps=convert(Int, ceil(abs( final_time / log_dt)))


root = div(N,2)
if MPI.Comm_rank(comm) == root
    print(" Running on $(MPI.Comm_size(comm)) processes\n")
end
rank=MPI.Comm_rank(comm)
MPIsize=MPI.Comm_size(comm)

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
println("got gates")

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

function propagate!(ψ,PH,propfunc)
  #ldims=linkdims(ψ)
  ψ,PH=propfunc(ψ,PH)
  return ψ,PH
end

opsum=tfd_holstein(N;omega=1.0,t=1.0,alpha=1.0, T=0.4, order=["phys_bos","el","anc_bos"])
states = [n == rank+1 ? [1,"Occ",1] : [1,"Emp",1] for n=1:N]
states=reduce(vcat,states)
@show states
@show typeof(sites[1])
ψ = MPS(ComplexF64,sites, states)
Nphys_sites=div(length(ψ),3)
H = MPO(opsum,sites)


PH=ProjMPO(H)
ψ0=deepcopy(ψ)
ϕ=nothing
##compute expectation values at t=0
if rank==root
    println("t=0 observables")
end
#ϕ0 = nothing
#ϕ0 = MPI.bcast(ψ0,root,comm)
ϕ = MPI.bcast(ψ,root,comm)
gf = measure_corr(ψ,ϕ,ψ0)
E = inner(ψ',H,ψ)
allgf=MPI.Gather(gf, root,comm)
allEs=MPI.Gather(E, root,comm)
if rank==root
    log(outf,"E",0.0,allEs)
    log(outf,"gf",0.0,allgf)
end
if rank==root
    println("starting prop")
end
@time begin
for i in range(1,total_Nsteps)
  if rank==root
    ψ,PH=propagate!(ψ,PH,tebd_step!)
  else
    ψ,PH=propagate!(ψ,PH,tebd_step!)
  end
  allgf,allgf2=measure_corr(ψ,ϕ,ψ0)
  @show allgf

  #ϕ = MPI.bcast(ψ,root,comm)
  #
  #gf = measure_corr(ψ,ϕ,conjugate=true)
  #gf2=measure_corr(ψ0,ϕ,conjugate=true)
  E = inner(ψ',H,ψ)
  #allgf=MPI.Gather(gf, root,comm)
  #allgf2=MPI.Gather(gf2, root,comm)
  
  #allEs=MPI.Gather(E, root,comm)
  if rank==root
    log(outf,"E",log_dt*i,allEs)
    log(outf,"gf",log_dt*i*2,allgf)
    log(outf,"gf2",log_dt*i,allgf2)
    
  end
end
end
#PH=ProjMPO(H)
end


  
