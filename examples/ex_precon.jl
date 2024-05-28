
using AtomsBase, DecoratedParticles, AtomsBuilder, 
      GeomOpt, Test, StaticArrays, Unitful, LinearAlgebra , 
      EmpiricalPotentials, SparseArrays, SaddleSearch

using AtomsCalculators: virial, forces, potential_energy
GO = GeomOpt      
DP = DecoratedParticles

import Random
Random.seed!(100)


##
 
sys = rattle!(bulk(:Si, cubic=true) * (20,2,2), 0.1)
sys = deleteat!(sys, 1)  # create a defect? 
@show length(sys)

sys0 = sys # AosSystem(sys)
sw = StillingerWeber()

# 1. Optimize without preconditioner (or, P = I)
sys1, result = GO.minimise(sys0, sw; g_tol = 1e-5)
# 219 force calls

# 2. Optimize with static preconditioner
P = GO.preconditioner(sys0, sw)
sysPst, resultPst = GO.minimise(sys0, sw; g_tol = 1e-5, precond = P)
# 21 force calls 

# 3. Optimize with dynamic preconditioner
#    this is much more expensize, works well with significant changes in the system 
sysPdyn, resultPdyn = GO.minimise(sys0, sw; g_tol = 1e-5, 
                        precond = sys -> GO.preconditioner(sys, sw))

##                                   
@info("     Without P : # iterations = $(result.iterations)")
@info(" With static P : # iterations = $(resultPst.iterations)")
@info("With dynamic P : # iterations = $(resultPdyn.iterations)")
