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
    
  function U_elel(siteinds, i::Int,j::Int; τ,alpha,omega,t,T) ##this is not parallelizable, would need to further split by even/odd terms
        elpos=2
        ancpos=3
        physpos=1
        
        s1=siteinds[3*(i-1)+elpos]
        s2=siteinds[3*(j-1)+elpos]
        h =
          (-t) * op("Cdag", s1)* op("C", s2) +
          (-t) * op("Cdag", s2)* op("C", s1)
          return exp(τ * h)
        end


function _get_gates_2ndorder(sites, tau,alpha,omega,t,T)
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
  end

function tfd_holstein_gates(sites, tau,order, alpha,omega,t,T)  ###electron hopping is not parallelizable in this form
    if order==2
      gates = _get_gates_2ndorder(sites, tau,alpha,omega,t,T)
      gates = collect(Iterators.flatten([gates,]))
      @show typeof(gates)
    elseif order==4 ###does not seem to perform as well as expected
      
      variant=1
      if variant==1
        s = 1.0 / ( 2.0 - 2.0^(1.0/3.0) )
        A = _get_gates_2ndorder(sites, s*tau,alpha,omega,t,T)
        B = _get_gates_2ndorder(sites, (1.0 -2.0*s)*tau,alpha,omega,t,T)
        gates = collect(Iterators.flatten([A,B,A]))
      elseif variant==2
        s2 = 1.0 / ( 4.0 - 4.0^(1.0/3.0) )
        A = _get_gates_2ndorder(sites, s2*tau,alpha,omega,t,T)
        B = _get_gates_2ndorder(sites, (1.0 - 4.0*s2)*tau,alpha,omega,t,T)
        gates = collect(Iterators.flatten([A,A,B,A,A]))
      end
      @show typeof(gates)
    end
    return gates
    
end