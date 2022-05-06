using ITensors
using ITensorTDVP
import MPI

MPI.Init()
comm = MPI.COMM_WORLD
N = 5
root = 0

if MPI.Comm_rank(comm) == root
    print(" Running on $(MPI.Comm_size(comm)) processes\n")
end
MPI.Barrier(comm)

if MPI.Comm_rank(comm) == root    
    N = 20
    s = siteinds(2,N)
    chi = 64
    psi = randomMPS(s;linkdims=chi)
else
    psi = nothing
end
let
i=MPI.Comm_rank(comm)
allranks = MPI.Gather(i,root,comm)
@show allranks
if MPI.Comm_rank(comm) == root
    @show allranks
end
for i in range(1,2)
    @time phi = MPI.bcast(psi, root, comm)
end
#println(maxlinkdim(phi)," ",MPI.Comm_rank(comm))
MPI.Finalize()
end