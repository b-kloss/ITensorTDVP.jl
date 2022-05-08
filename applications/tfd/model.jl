#include("boson.jl")#
function theta(omega, beta)
    return atanh(exp(-beta * omega / 2.))
end

function V(alpha,omega,T)
    if T==0.0
        return alpha
    else
        return alpha*cosh(theta(omega,1.0/T))
    end
end

function Vtilde(alpha,omega,T)
    if T==0.0
        return 0.0
    else
        return alpha*sinh(theta(omega,1.0/T))
    end
end
    

function tfd_holstein(n;omega=1.0,t=1.0, alpha=1.0, T=0.0, order=["phys_bos","el","anc_bos"])
    os = OpSum()
    @assert order==["phys_bos","el","anc_bos"]
    elpos=2
    ancpos=3
    physpos=1
    for j in 1:(n - 1)
      elj=elpos
      
      os += -t,"Cdag", 3*(j-1)+elj, "C", 3*(j) +elj
      os += -t,"Cdag", 3*(j) +elj , "C", 3*(j-1)+elj
    end
    for j in 1:n
      os += Vtilde(alpha,omega,T),"n", 3*(j-1)+elpos, "A", 3*(j-1) + ancpos   #Vdag
      os += Vtilde(alpha,omega,T),"n", 3*(j-1)+elpos, "Adag", 3*(j-1) + ancpos
      os += V(alpha,omega,T),"A", 3*(j-1) + physpos, "n", 3*(j-1)+elpos   #V
      os += V(alpha,omega,T), "Adag", 3*(j-1) + physpos, "n", 3*(j-1)+elpos
      os += omega, "N", 3*(j-1)+physpos #local oscillator
      os += -omega, "N", 3*(j-1)+ancpos #local ancilla oscillator
    end
    return os
  end

