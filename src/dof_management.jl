
using StaticArrays, AtomsBase, DecoratedParticles, Unitful 
using LinearAlgebra: I 

"""
`DofManager`: 

Constructor:
```julia
DofManager(sys; free=..., clamp=..., mask=...)
```
Set at most one of the kwargs:
* no kwarg: all atoms are free
* `free` : list of free atom indices (not dof indices)
* `clamp` : list of clamped atom indices (not dof indices)
* `mask` : 3 x N Bool array to specify individual coordinates to be clamped

### Meaning of dofs

On call to the constructor, `DofManager` stores positions and cell
`X0, C0`, dofs are understood *relative* to this initial configuration.
`dofs(sys, dm::DofManager)` returns a vector that represents the displacement 
and a deformation matrix `(U, F)`. The new configuration extracted from a dof vector 
is understood as 
* The new cell: `C = F * C0` 
* The new positions: `𝐫[i] = F * (X0[i] + U[i])`
One aspect of this definition is that clamped atom positions still change via
the deformation `F`. This is natural in the context of optimizing the 
cell shape. 
"""
mutable struct DofManager{D,T}
   variablecell::Bool
   ifree::Vector{Int}   # extract the free position dofs 
   X0::Vector{SVector{D,T}}       # reference positions
   C0::NTuple{D, SVector{D, T}}   # reference cell 
end


# NOTES: 
#  - the reference positions maybe are not needed 
#       and could be removed in the next iteration
#  - length units are implicitly given by the units in X0 and C0. 
#    unit stripping and conversions will be careful to maintin consistency  

# ========================================================================
#  Constructors 

function DofManager(sys::AbstractSystem{D};  
                    variablecell = false, 
                    free = nothing, 
                    clamp = nothing, 
                    mask = nothing,  )  where {D} 
   if D != 3 
      error("this package assumes d = 3; please file an issue if you neeed a different use case")                    
   end 
   X0 = position(sys)
   C0 = tuple(bounding_box(sys)...)   
   ifree = analyze_mask(sys, free, clamp, mask)
   return DofManager(variablecell, ifree, X0, C0)
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

length_unit(dm::DofManager) = unit(dm.X0[1][1])
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


function _pos2dofs(U, dofmgr) 
   TFL = eltype(ustrip(U[1]))
   return reinterpret(TFL, U)[dofmgr.ifree]
end

function _dofs2pos(x::AbstractVector{T}, dofmgr)  where {T} 
   TU = eltype(dofmgr.X0) 
   u = zeros(T, 3 * length(dofmgr.X0))
   @show length(dofmgr.ifree) 
   @show length( _posdofs(x, dofmgr))
   u[dofmgr.ifree] .= _posdofs(x, dofmgr)
   return reinterpret(TU, u)
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
   TFL = eltype(ustrip(X[1]))

   if fixedcell(dofmgr)
      # there is an allocation here that could maybe be avoided 
      return _pos2dofs(X - dofmgr.X0, dofmgr)::Vector{TFL}
   end

   # variable cell case: note we already checked units and can strip 
   bb = bounding_box(sys) 
   F = ustrip.(hcat(bb...)) / ustrip.(hcat(dofmgr.C0...))
   # Xi = F * (X0i + Ui)  =>  Ui = F \ Xi - X0i
   U = [ F \ X[i] - dofmgr.X0[i] for i = 1:length(X) ]
   return [ _pos2dofs(U, dofmgr); 
            _defm2dofs(F, dofmgr) ]::Vector{TFL} 
end


function set_dofs!(sys::AbstractSystem, dofmgr::DofManager, 
                   x::AbstractVector{T} ) where {T <: AbstractFloat}
   check_length_units(sys, dofmgr)

   # get the displacement from the dof vector    
   U, F = _dofs2posdefm(x, dofmgr)

   # convert the displacements to positions 
   X = [ F * (dofmgr.X0[i] + U[i]) for i = 1:length(U) ]
   bb_old = bounding_box(sys)
   bb_new = ntuple(i -> F * bb_old[i], 3)
   # and update the system 
   set_positions!(sys, X)
   set_bounding_box!(sys, bb_new)
   return sys
end



# ========================================================================
#   Compute the gradient with respect to dofs 
#   from forces and virials 


function energy_dofs(sys, calc, dofmgr, x)

end

function gradient_dofs(sys, calc, dofmgr, x)

end


# function gradient(calc, at::Atoms)
#    if fixedcell(at)
#       return rmul!(mat(forces(calc, at))[at.dofmgr.xfree], -1.0)
#    end
#    F = cell(at)'
#    A = F * inv(at.dofmgr.F0)
#    G = forces(calc, at)
#    for n = 1:length(G)
#       G[n] = - A' * G[n]
#    end
#    S = - virial(calc, at) * inv(F)'                  # ∂E / ∂F
#    # S += at.dofmgr.pressure * sigvol_d(at)'     # applied stress  TODO: revive this!
#    return [ mat(G)[at.dofmgr.xfree]; Array(S)[:] ]
# end

