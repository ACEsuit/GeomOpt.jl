#
# WARNING : this example runs but gives 
#           unphysical results
#

using AtomsBase, DecoratedParticles, AtomsBuilder, 
      GeomOpt, Test, StaticArrays, Unitful, LinearAlgebra , 
      EmpiricalPotentials, SparseArrays, SaddleSearch

using AtomsCalculators: virial, forces, potential_energy
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

# pick an arbitrary FCC system (needs to be closed-packed, hcp would also be ok)
r0 = AtomsBuilder.Chemistry.rnn(:Si)
sys0 = AosSystem( bulk(:Si, cubic=true) * 3 )
# x2 is a nearest-neighbour of x1
@assert norm(position(sys0, 1)) == 0u"Å"
@assert norm(position(sys0, 2)) <= 1.0001 * r0
𝐫2 = position(sys0, 2)
# delete the second atom, then make a copy of the system 
# and move x1 into the position of what was previously x2.
deleteat!(sys0.particles, 2)
sys0 = AosSystem(sys0)
calc = StillingerWeber()

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

## run a NEB method 
preconI = SaddleSearch.localPrecon(precon = [I], precon_prep! = (P, x) -> P)
# neb = ODENEB(reltol=1e-2, k=0.0002, interp=3, tol = 1e-2, maxnit = 100,
#               precon_scheme = preconI, verbose = 1)
# xx_neb, log_neb, _ = run!(neb, obj_f, obj_g, Path(xx_init))
# @show log_neb[:maxres][end]

stringmethod = ODEString(reltol=0.1, tol = 1e-2, maxnit = 100, 
                         precon_scheme = preconI, verbose = 1)
xx_string, log_string, _ = SaddleSearch.run!(stringmethod, obj_f, obj_g, Path(xx_init))
@show log_string[:maxres][end]


# same but now with preconditioning  (still doesn't converge ...)
function obj_precon(x)
   GeomOpt.set_dofs!(sys0, dofmgr, x)
   preconditioner(sys0, calc)
end
preconP = SaddleSearch.localPrecon(
               precon = obj_precon.(xx_init), 
               precon_prep! = (P, xx) -> obj_precon.(xx) )

# neb_P = ODENEB(reltol=1e-2, k=0.0002, interp=3, tol = 1e-2, maxnit = 100,
#               precon_scheme = preconP, verbose = 1)
# xx_neb_P, log_neb_P, _ = run!(neb_P, obj_f, obj_g, Path(xx_init))
# @show log_neb_P[:maxres][end]

string_P = ODEString(reltol=0.1, tol = 1e-2, maxnit = 100, 
                         precon_scheme = preconP, verbose = 1)
xx_string_P, log_string_P, _ = SaddleSearch.run!(stringmethod, obj_f, obj_g, Path(xx_init))
@show log_string_P[:maxres][end]

##

@info("energy along initial and neb paths")
# E_neb = [ obj_f(x) for x in xx_neb ]
# δE_neb = round.(E_neb .- E_neb[1], digits=3)*u"eV"
# E_neb_P = [ obj_f(x) for x in xx_neb_P ]
# δE_neb_P = round.(E_neb_P .- E_neb_P[1], digits=3)*u"eV"
# display( hcat(δE_init, δE_neb, δE_neb_P) )

E_string = [ obj_f(x) for x in xx_string ]
δE_string = round.(E_string .- E_string[1], digits=3)*u"eV"
E_string_P = [ obj_f(x) for x in xx_string_P ]
δE_string_P = round.(E_string_P .- E_string_P[1], digits=3)*u"eV"
display( hcat(δE_init, δE_string, δE_string_P) )


