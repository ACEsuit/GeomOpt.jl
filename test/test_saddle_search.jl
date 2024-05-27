
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


#

# first get the ground state unit cell 
# this needs some fixes in DP
# ucell0 = AosSystem(bulk(:Si, cubic=true))
# sw = StillingerWeber() 
# ucell1, _ = GO.minimise(ucell0, sw; g_tol = 1e-5, variablecell=true)

r0 = AtomsBuilder.Chemistry.rnn(:Cu)
sys0 = AosSystem( bulk(:Cu, cubic=true) * 3 )
# x2 is a nearest-neighbour of x1
@assert norm(position(sys0, 1)) == 0u"Å"
@assert norm(position(sys0, 2)) <= 1.0001 * r0
𝐫2 = position(sys0, 2)
# delete the second atom and move at1 onto the midpoint 
# the new system sys1 is an initial guess for a saddle point 
deleteat!(sys0.particles, 2)
sys1 = SoaSystem(sys0)
set_position!(sys1, 1, 0.5 * 𝐫2)

# LJ is not a good model for Cu, but a perfectly ok model for a generic 
# FCC crystal at 0T 
zU = sys0[1].𝑍
calc = LennardJones(-1.0u"eV", r0, zU, zU, 3*r0,)

# create a system -> dof mapping 
dofmgr1 = GeomOpt.DofManager(sys1; variablecell=false)
obj_f, obj_g = GeomOpt.get_obj_fg(sys1, calc, dofmgr1)
x1 = GeomOpt.get_dofs(sys1, dofmgr1)

# to get the initial orientation of the dimer, we get the reference position 
# first and then take the difference 
x0 = GeomOpt.get_dofs(sys0, dofmgr1)
v0 = x1 - x0  # initial rotation 

dimer = StaticDimer(a_trans=0.001, a_rot=0.001, len=1e-3, maxnumdE=3000)
x, v, report = SaddleSearch.run!(dimer, obj_f, obj_g, x1, v0)

obj_f(x1)
norm(obj_g(x1))
