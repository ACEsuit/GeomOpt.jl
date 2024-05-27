
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

# LJ is not a good model for Cu, but a perfectly ok model for a generic 
# FCC crystal at 0T 
zU = atomic_number(sys0, 1)
calc = LennardJones(-1.0u"eV", r0, zU, zU, 3*r0,)

# generate an initial guess for the MEP 
Nimg = 11
path_init = [ (sys = deepcopy(sys0); set_position!(sys, 1, t * 𝐫2); sys)
               for t in range(0.0, 1.0, length = Nimg) ]
E_init = potential_energy.(path_init, Ref(calc))
@info("Initial guess, energy difference")
display(round.(u"eV", E_init .- E_init[1], digits=2))

##
# create a system -> dof mapping 
dofmgr = GeomOpt.DofManager(sys0; variablecell=false)
obj_f, obj_g = GeomOpt.get_obj_fg(sys0, calc, dofmgr)
xx_init = [ GeomOpt.get_dofs(sys, dofmgr) for sys in path_init ]

# just a quick sanity check - note this will fail if we introduce 
# more general energy-nondimensionalization procedures 
E_check = [ obj_f(x) for x in xx_init ]
all(E_check .≈ ustrip.(E_init))

## run a String method 
# this doesn't converge, too few iterations, but demonstrates usage

preconI = SaddleSearch.localPrecon(precon = [I], precon_prep! = (P, x) -> P)
stringmethod = ODEString(reltol=0.1, tol = 1e-2, maxnit = 100, 
                         precon_scheme = preconI, verbose = 1)
PATHx, PATHlog, _ = SaddleSearch.run!(stringmethod, obj_f, obj_g, Path(xx_init))
@show PATHlog[:maxres][end]

@info("energy along initial / final path")
E_final = [ obj_f(x) for x in PATHx ]
display( hcat(round.(u"eV", E_init .- E_init[1], digits=2), 
              round.(E_final .- E_final[1], digits=2)*u"eV") )

              