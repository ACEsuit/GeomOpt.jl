
using StaticArrays, AtomsBase, DecoratedParticles, Unitful 
using LinearAlgebra: I, mul! 

"""
`DofManager`: 

Constructor:
```julia
DofManager(sys; variablecell = false, r0 =..., free=..., clamp=..., mask=...)
```
* `variablecell` determines whether the cell is fixed or allowed to change 
  during optimization 
* `r0` is a reference length-scale, default set to one (in the unit of sys), 
   this is used to non-dimensionalize the degrees of freedom. 

In addition set at most one of the kwargs:
* no kwarg: all atoms are free
* `free` : list of free atom indices (not dof indices)
* `clamp` : list of clamped atom indices (not dof indices)
* `mask` : 3 x N Bool array to specify individual coordinates to be clamped

### Meaning of dofs

On call to the constructor, `DofManager` stores positions and cell
`X0, C0`, dofs are understood *relative* to this initial configuration.
`get_dofs(sys, dm::DofManager)` returns a vector that represents the
non-dimensional displacement and a deformation matrix `(U, F)`. The new configuration extracted from a dof vector 
is understood as 
* The new cell: `C = F * C0` 
* The new positions: `𝐫[i] = F * (X0[i] + U[i] * r0)`
One aspect of this definition is that clamped atom positions still change via
the deformation `F`. This is natural in the context of optimizing the 
cell shape. 
"""
mutable struct DofManager{D,T}
   variablecell::Bool
   ifree::Vector{Int}   # extract the free position dofs 
   r0::T 
   X0::Vector{SVector{D,T}}       # reference positions
   C0::NTuple{D, SVector{D, T}}   # reference cell 
end


# NOTES: 
#  - length units are implicitly given by the units in X0, C0, r0. 
#    there should be no explicit length-unit stripping but this should be 
#    implicity through the reference length-scale r0
#  - at the moment energy-nondimensionalization is achieved simply by 
#    stripping. A better approach would be to enforce this to happen in 
#    the preconditioner, which could simply be a rescaling operation. 

# ========================================================================
#  Constructors 

function DofManager(sys::AbstractSystem{D};  
                    variablecell = false, 
                    r0 = _auto_r0(sys), 
                    free = nothing, 
                    clamp = nothing, 
                    mask = nothing,  )  where {D} 
   if D != 3 
      error("this package assumes d = 3; please file an issue if you neeed a different use case")                    
   end 
   X0 = position(sys)
   C0 = tuple(bounding_box(sys)...)   
   ifree = analyze_mask(sys, free, clamp, mask)
   return DofManager(variablecell, ifree, r0, X0, C0)
end

function _auto_r0(sys)
   r = position(sys, 1)[1] 
   return r0 = one(ustrip(r)) * unit(r) 
end 

"""
`analyze_mask` : helper function to generate list of dof indices from
lists of atom indices indicating free and clamped atoms
"""
function analyze_mask(sys, free, clamp, mask)
   if length(findall((free != nothing, clamp != nothing, mask != nothing))) > 1
      error("DofManager: only one of `free`, `clamp`, `mask` may be provided")
   elseif all( (free == nothing, clamp == nothing, mask == nothing) )
      # in this case (default) all atoms are free
      return collect(1:3*length(sys))
   end
   # determine free dof indices
   Nat = length(sys)
   if clamp != nothing
      # revert to setting free
      free = setdiff(1:Nat, clamp)
   end
   if free != nothing
      # revert to setting mask
      mask = fill(false, 3, Nat)
      if !isempty(free)
         mask[:, free] .= true
      end
   end
   return findall(mask[:])
end


# ========================================================================
#   DOF Conversions

length_unit(dm::DofManager) = unit(dm.r0)
length_unit(sys::AbstractSystem) = unit(position(sys, 1)[1])

function check_length_units(sys, dm::DofManager) 
   if length_unit(dm) != length_unit(sys) 
      error("System `sys` and DofManager have inconsistent units.")
   end 
   if length(sys) != length(dm.X0)
      error("System `sys` and DofManager have inconsistent size.")
   end
   nothing 
end 

fixedcell(dofmgr::DofManager) = 
      !dofmgr.variablecell 

variablecell(dofmgr::DofManager) = 
      dofmgr.variablecell 

# there is a type-instability here!! 
_posdofs(x, dofmgr::DofManager) = 
      dofmgr.variablecell ? (@view x[1:end-9]) : x 


_pos2dofs(U::AbstractVector{SVector{3, T}}, dofmgr)  where {T} = 
      @view(reinterpret(T, U)[dofmgr.ifree])

function _dofs2pos(x::AbstractVector{T}, dofmgr)  where {T} 
   u = zeros(T, 3 * length(dofmgr.X0))
   u[dofmgr.ifree] .= _posdofs(x, dofmgr)
   return reinterpret(SVector{3, T}, u)
end

_celldofs(x, dofmgr) = 
      dofmgr.variablecell ? x[end-8:end] : missing 

function _defm2dofs(F, dofmgr)
   return Matrix(F)[:]
end

# there is another type-instability here!! 
function _dofs2defm(x::AbstractVector{T}, dofmgr) where {T}
   if dofmgr.variablecell 
      return SMatrix{3, 3}(_celldofs(x, dofmgr))
   else 
      return SMatrix{3, 3}( T[1 0 0; 0 1 0; 0 0 1] )
   end 
end

_dofs2posdefm(x, dofmgr) = 
      _dofs2pos(x, dofmgr), _dofs2defm(x, dofmgr)


function get_dofs(sys::AbstractSystem, dofmgr::DofManager)
   check_length_units(sys, dofmgr)

   # obtain the positions and their underlying floating point type 
   X = position(sys)

   if fixedcell(dofmgr)
      # there are allocations here that could maybe be avoided 
      return collect(_pos2dofs((X - dofmgr.X0)/dofmgr.r0, dofmgr))
   end

   # variable cell case: note we already checked units and can strip 
   #   (otherwise we would have problems inverting)
   bb = bounding_box(sys) 
   F = ustrip.(hcat(bb...)) / ustrip.(hcat(dofmgr.C0...))
   # Xi = F * (X0i + Ui * r0)  =>  Ui = (F \ Xi - X0i) / r0 
   U = [ (F \ X[i] - dofmgr.X0[i]) / dofmgr.r0 for i = 1:length(X) ]
   return [ _pos2dofs(U, dofmgr); 
            _defm2dofs(F, dofmgr) ]
end


function set_dofs!(sys::AbstractSystem, dofmgr::DofManager, 
                   x::AbstractVector{T} ) where {T <: AbstractFloat}
   check_length_units(sys, dofmgr)

   # get the displacement from the dof vector    
   U, F = _dofs2posdefm(x, dofmgr)

   # convert the displacements to positions 
   X = [ F * (dofmgr.X0[i] + U[i] * dofmgr.r0) for i = 1:length(U) ]
   bb_old = dofmgr.C0
   bb_new = ntuple(i -> F * bb_old[i], 3)
   # and update the system 
   set_bounding_box!(sys, bb_new)
   set_positions!(sys, X)
   return sys
end



# ========================================================================
#   Compute the gradient with respect to dofs 
#   from forces and virials 

using AtomsCalculators: potential_energy, energy_forces_virial, 
            forces, virial 

function energy_dofs(sys, calc, dofmgr, x)
   set_dofs!(sys, dofmgr, x)
   return ustrip(potential_energy(sys, calc))
end


function gradient_dofs(sys, calc, dofmgr, x::AbstractVector{T}) where {T} 
   set_dofs!(sys, dofmgr, x)

   EFV = energy_forces_virial(sys, calc) 
   frc = EFV.forces 
   vir = EFV.virial 

   # now transform forces and virial into a gradient w.r.t. x 
   # ------------------------------------------------------------
   # fixed cell version 
   # fi = - ∇_𝐫i E  [eV/A]
   # 𝐫i = X0[i] + r0 * U[i] 
   # g_iα = - fiα * r0  [eV] => same unit as E so can strip 
   if fixedcell(dofmgr)
      g = [ ustrip( - dofmgr.r0 * f ) for f in frc ]
      return collect(_pos2dofs(g, dofmgr))::Vector{T}
   end

   # variable cell version 
   # fi = - ∇_𝐫i E  [eV/A]     𝐫i = F * (X0[i] + r0 * U[i])
   # ∇_𝐮i' = - fi' * ∂𝐫i/∂𝐮i = - fi' * (r0 * F)   =>   ∇_𝐮i = - F' * r0 * fi 
   # ∂F E = - virial
   F = _dofs2defm(x, dofmgr)
   g_pos = [ - ustrip(dofmgr.r0 * F' * f) for f in frc ]

   return [ _pos2dofs(g_pos, dofmgr); 
            ( - ustrip.(vir) / F' )[:] ]::Vector{T}
end

