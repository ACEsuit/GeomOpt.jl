
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
δE_init = round.(u"eV", E_init .- E_init[1], digits=2)
display(δE_init)

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
xx_string, log_string, _ = SaddleSearch.run!(stringmethod, obj_f, obj_g, Path(xx_init))
@show log_string[:maxres][end]

@info("energy along initial / final path")
E_string = [ obj_f(x) for x in xx_string ]
δE_string = round.(E_string .- E_string[1], digits=2)*u"eV"
display( hcat(δE_init, δE_string) )

## Now try a NEB               

neb = ODENEB(reltol=1e-2, k=0.0002, interp=3, tol = 1e-2, maxnit = 100,
              precon_scheme = preconI, verbose = 1)
xx_neb, log_neb, _ = run!(neb, obj_f, obj_g, Path(xx_init))
@show log_neb[:maxres][end]
    
@info("energy along initial / string / neb paths")
E_neb = [ obj_f(x) for x in xx_neb ]
δE_neb = round.(E_neb .- E_neb[1], digits=2)*u"eV"
display( hcat(δE_init, δE_string, δE_neb) )

