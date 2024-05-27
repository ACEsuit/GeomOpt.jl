
using AtomsBase, DecoratedParticles, AtomsBuilder, 
      GeomOpt, Test, StaticArrays, Unitful, LinearAlgebra , 
      EmpiricalPotentials, SaddleSearch

using AtomsCalculators: virial, forces
GO = GeomOpt      
DP = DecoratedParticles

import Random
Random.seed!(100)

##

# first get the ground state unit cell 
# this needs some fixes in DP
# ucell0 = AosSystem(bulk(:Si, cubic=true))
# sw = StillingerWeber() 
# ucell1, _ = GO.minimise(ucell0, sw; g_tol = 1e-5, variablecell=true)

sys0 = AosSystem( bulk(:Si, cubic=true) * 2 )
# sys0[2] is a neighbour of sys0[1] 
@assert norm(position(sys0, 1)) == 0u"Å"
@assert norm(position(sys0, 2)) < 2.36u"Å"
𝐫2 = position(sys0, 2)
# delete the second atom and move at1 onto the midpoint 
deleteat!(sys0.particles, 2)
sys1 = SoaSystem(sys0)
set_position!(sys1, 1, 0.5 * 𝐫2)

# create a system -> dof mapping 
dofmgr1 = GeomOpt.DofManager(sys1; variablecell=false)
obj_f, obj_g = GeomOpt.get_obj_fg(sys1, sw, dofmgr1)
x1 = GeomOpt.get_dofs(sys1, dofmgr1)

x0 = GeomOpt.get_dofs(sys0, dofmgr)
v0 = x1 - x0  # initial rotation 

dimer = StaticDimer(a_trans=0.0001, a_rot=0.0001, len=1e-3, maxnumdE=3000)
x, v, report = SaddleSearch.run!(dimer, obj_f, obj_g, x1, v0)
