


using AtomsBase, DecoratedParticles, AtomsBuilder, 
      GeomOpt, Test, StaticArrays, Unitful, LinearAlgebra , 
      EmpiricalPotentials

using AtomsCalculators: virial, forces
GO = GeomOpt      
DP = DecoratedParticles

import Random
Random.seed!(100)

##
#
# Test 1 : equilibrate a perturbed Si crystal with frozen cell 
#          the positions should move back into crystalline 
#          hence up to a shift they must be the same as sys2
#
sys0 = AosSystem( bulk(:Si, cubic=true) )
sys1 = AosSystem( rattle!(bulk(:Si, cubic=true), 0.1) )
sw = StillingerWeber()
sys2, result = GO.minimise(sys1, sw; g_tol = 1e-5)

@test result.g_residual < 1e-5
@test bounding_box(sys2) == bounding_box(sys1) == bounding_box(sys0)
X0 = [ x - position(sys0, 1) for x in position(sys0, :) ]
X2 = [ x - position(sys2, 1) for x in position(sys2, :) ]
@test maximum(ustrip.(norm.(X2 .- X0)) .< 1e-6)

##
#
# Test 2 : equilibrate a unitcell 
#
sys1 = AosSystem( bulk(:Si) )
sw = StillingerWeber()
sys2, result = GO.minimise(sys1, sw; variablecell=true, 
                                     g_tol = 1e-5)
@test result.g_residual < 1e-5
@test norm( ustrip.(virial(sys2, sw)) ) < 1e-5

##
#
# Test 3 : equilibrate positions and cell
#
sys1 = AosSystem( rattle!(bulk(:Si, cubic=true)*2, 0.1) )
sw = StillingerWeber()
sys2, result = GO.minimise(sys1, sw; variablecell=true, 
                                     g_tol = 1e-3)
@test result.g_residual < 1e-3
f = ustrip.(forces(sys2, sw))
@test maximum(norm.(f)) < 1e-3
@test norm( ustrip.(virial(sys2, sw)) ) < 2e-3
# 2e-3 is ok here because the virial is not actually the gradient 
# the relationships is a little tricky

