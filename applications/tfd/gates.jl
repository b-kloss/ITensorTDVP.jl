op(s::String,i::Index) = ITensors.op(s::String,i::Index)

function U_elph(siteinds, i::Int; τ,alpha,omega,t,T)
    elpos=2
    ancpos=3
    physpos=1
    
    s1=siteinds[3*(i-1)+elpos]
    s2=siteinds[3*(i-1)+physpos]
    h =
        V(alpha,omega,T) * op("n", s1) * op("A",s2) +
        V(alpha,omega,T) * op("n", s1) * op("Adag", s2) +
        omega * op("Id",s1) * op("N",s2)
      return exp(τ * h)
    end
  
  function U_elanc(siteinds, i::Int64; τ,alpha,omega,t,T)
      elpos=2
      ancpos=3
      physpos=1
      
      s1=siteinds[3*(i-1)+elpos]
      s2=siteinds[3*(i-1)+ancpos]
      h =
          Vtilde(alpha,omega,T)*op("n", s1) * op("A",s2) +
          Vtilde(alpha,omega,T)*op("n", s1) * op("Adag", s2) +
          (-omega) * op("Id",s1) * op("N",s2)
        return exp(τ * h)
      end
    
  function U_elel(siteinds, i::Int,j::Int; τ,alpha,omega,t,T)
        elpos=2
        ancpos=3
        physpos=1
        
        s1=siteinds[3*(i-1)+elpos]
        s2=siteinds[3*(j-1)+elpos]
        h =
          t* op("Cdag", s1)* op("C", s2) +
          t* op("Cdag", s2)* op("C", s1)
          return exp(τ * h)
        end

function tfd_holstein_gates(sites, tau,order, alpha,omega,t,T)
    N=div(length(sites),3)
    anc_gates = Vector{ITensor}()
    phys_gates = Vector{ITensor}()
    el_gates = Vector{ITensor}()
    for i in 1:N
        push!(anc_gates,U_elanc(sites,i;τ=-tau * im /2., alpha,omega,t,T))
        push!(phys_gates,U_elph(sites,i;τ=-tau * im /2., alpha,omega,t,T))
    end
    for i in 1:N-1
        push!(el_gates,U_elel(sites,i,i+1;τ=-tau * im /2., alpha,omega,t,T))
    end
    gates=[phys_gates,anc_gates,el_gates,el_gates,reverse(anc_gates),reverse(phys_gates)]
    return gates
    #anc_gates=ops([("expτ-elanc", (n), (τ=-tau * im /2., alpha=alpha,omega=omega,t=t,T=T)) for n in range(1,N,step=1)], sites)
     
end