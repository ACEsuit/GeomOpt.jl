
using EmpiricalPotentials, SparseArrays

function preconditioner(sys::AbstractSystem, calc)
   Nat = length(sys)
   P = zeros(3*Nat, 3*Nat) 
   nlist = EmpiricalPotentials.PairList(sys, cutoff_radius(calc))
   for i = 1:Nat 
      Js, Rs, Zs, z0 = EmpiricalPotentials.get_neighbours(sys, calc, nlist, i)
      Pi = EmpiricalPotentials.precon(calc, Rs, Zs, z0) 

      Nr = length(Js)
      D = 3 
      Ji = (i - 1) * D .+ (1:D)
      for (α1, j1) in enumerate(Js), (α2, j2) in enumerate(Js)
         # A1 = (α1-1) * D .+ (1:D)
         # A2 = (α2-1) * D .+ (1:D)
         J1 = (j1-1) * D .+ (1:D)
         J2 = (j2-1) * D .+ (1:D)
         P[J1, J2] += Pi[α1, α2]
         P[J1, Ji] -= Pi[α1, α2]
         P[Ji, J2] -= Pi[α1, α2]
         P[Ji, Ji] += Pi[α1, α2]
      end
   end
   return sparse(P) + 0.01*I 
end
