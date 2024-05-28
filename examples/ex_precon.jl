using AtomsBase, DecoratedParticles, AtomsBuilder, 
      GeomOpt, Test, StaticArrays, Unitful, LinearAlgebra , 
      EmpiricalPotentials, SparseArrays

using AtomsCalculators: virial, forces
GO = GeomOpt      
DP = DecoratedParticles

import Random
Random.seed!(100)

function preconditioner(sys::AbstractSystem, sw::StillingerWeber)
   Nat = length(sys)
   P = zeros(3*Nat, 3*Nat) 
   nlist = EmpiricalPotentials.PairList(sys, cutoff_radius(sw))
   for i = 1:Nat 
      Js, Rs, Zs, z0 = EmpiricalPotentials.get_neighbours(sys, sw, nlist, i)
      Pi = EmpiricalPotentials.precon(sw, Rs, Zs, z0) 

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

##

sys0 = AosSystem( rattle!(bulk(:Si, cubic=true) * (10,2,2), 0.1) )
@show length(sys)
sw = StillingerWeber()
sys1, result = GO.minimise(deepcopy(sys0), sw; g_tol = 1e-5)

P = preconditioner(sys0, sw)
sysP, resultP = GO.minimise(deepcopy(sys0), sw; g_tol = 1e-5, precond = P)

@info("Without P : # iterations = $(result.iterations)")
@info("   With P : # iterations = $(resultP.iterations)")

##

