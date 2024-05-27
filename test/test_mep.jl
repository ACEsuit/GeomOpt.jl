
using AtomsBase, DecoratedParticles, AtomsBuilder, 
      GeomOpt, Test, StaticArrays, Unitful, LinearAlgebra , 
      EmpiricalPotentials, SaddleSearch

using AtomsCalculators: virial, forces, potential_energy
GO = GeomOpt      
DP = DecoratedParticles

import Random
Random.seed!(100)

##

# pick an arbitrary FCC system (needs to be closed-packed, hcp would also be ok)
r0 = AtomsBuilder.Chemistry.rnn(:Cu)
sys0 = AosSystem( bulk(:Cu, cubic=true) * 3 )
# x2 is a nearest-neighbour of x1
@assert norm(position(sys0, 1)) == 0u"Å"
@assert norm(position(sys0, 2)) <= 1.0001 * r0
𝐫2 = position(sys0, 2)
# delete the second atom, then make a copy of the system 
# and move x1 into the position of what was previously x2.
deleteat!(sys0.particles, 2)
sys0 = SoaSystem(sys0)
sys1 = deepcopy(sys0)
set_position!(sys1, 1, 𝐫2)

# LJ is not a good model for Cu, but a perfectly ok model for a generic 
# FCC crystal at 0T 
zU = atomic_number(sys0, 1)
calc = LennardJones(-1.0u"eV", r0, zU, zU, 3*r0,)

# create a system -> dof mapping 
dofmgr = GeomOpt.DofManager(sys0; variablecell=false)
obj_f, obj_g = GeomOpt.get_obj_fg(sys0, calc, dofmgr1)
x0 = GeomOpt.get_dofs(sys0, dofmgr)
x1 = GeomOpt.get_dofs(sys1, dofmgr)

# x0 is just a bunch of zeros. (since it is still in the reference)
# but x1 will have the first three entries changed 
all(iszero, x0)         # true 
x1[1:3] == ustrip.(𝐫2)  # true 
norm(x1[4:end])         # 0.0 

# x0, x1 now represent the endpoints of the string. We now convert that 
# into a straight line segment in configuration space connecting those 
# two points 
Nimag = 20
sys_path_init = [ (sys = deepcopy(sys0); set_position!(sys, 1, t * 𝐫2); sys)
                  for t in range(0.0, 1.0, length = Nimag) ]
dof_path_init = [ GeomOpt.get_dofs(sys, dofmgr) 
                  for sys in sys_path_init ]         
# path_init = Path( [ (1-t) * x0 + t * x1
#                     for t in range(0, 1, length = Nimag) ] )

##
X = ustrip.(position(sys1))
rr = [ (i==j) * 100 + norm(X[i] - X[j]) for i = 1:length(X), j = 1:length(X)]
minimum(rr)

obj_f(x1)
potential_energy(sys1, calc)

# look at the energy along the initial path to check that we didn't do 
# anything stupid.                     
E_init = [ obj_f(x) for x in path_init.x  ]
@info("Energy along initial path")
display(round.(E_init, digits = 2))

# run a NEB and or 

path = ODEString(reltol=0.1, tol = 1e-2, maxnit = 1000, verbose = 1)
PATHx, PATHlog, _ = SaddleSearch.run!(path, obj_f, obj_g, path_init)
@test PATHlog[:maxres][end] <= path.tol



dimer = StaticDimer(a_trans=1e-3, a_rot=1e-3, len=1e-3, maxnumdE=1000, 
                    tol_trans = 1e-2, tol_rot = 1e-1)
x, v, report = SaddleSearch.run!(dimer, obj_f, obj_g, x1, v0)
@info("|∇E(x_final)|_inf = $(norm(obj_g(x), Inf))")

# TODO: add hessian spectrum test once relevant PRs are merged


bb = BBDimer(a0_trans=1e-3, a0_rot=1e-3, len=1e-3, maxnumdE=1000, 
                    tol_trans = 1e-2, tol_rot = 1e-1)
x, v, report = SaddleSearch.run!(bb, obj_f, obj_g, x1, v0)
@info("|∇E(x_final)|_inf = $(norm(obj_g(x), Inf))") 

# TODO: add hessian spectrum test once relevant PRs are merged