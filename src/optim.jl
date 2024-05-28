
import Optim
import LineSearches

using Optim: OnceDifferentiable, optimize, ConjugateGradient, LBFGS
using LineSearches: BackTracking

using LinearAlgebra: I

using AtomsBase, AtomsCalculators, DecoratedParticles


export minimise!


"""
`minimise!(sys, calc; kwargs...)`: minimize potential energy

## Keyword arguments:
* `precond = I` : preconditioner; more below
* `frc_tol = 1e-6u"eV/Å"` : force tolerance (max-norm)
* `vir_tol = 1e-6u"eV" : virial / Nat tolerance 
* `E_tol = 0.0u"eV"` : total energy tolerance`
* `Optimiser = :auto`, `:auto` should always work
* `verbose = 1`: 0 : no output, 1 : final, 2 : iteration and final
* `store_trace = false` : store history of energy and norm of forces
* `extended_trace = false`: also store full history of postions and forces
* `maxstep = Inf`: maximum step size, useful if initial gradient is very large
* `callback = nothing`: callback function to pass to `optimize()`, e.g. to use alternate convergence criteria
* `frc_call_limit = 1000` : set a limit on the number of force calls before giving up 

## Preconditioner

`precond` may be a valid preconditioner, e.g., `I` or `Exp(at)`, or one of
the following symbols

* `:auto` : the code will make the best choice it can with the available
   information
* `:exp` : will use `Exp(at)`
* `:id` : will use `I`
"""
function minimise(sys, calc;
                   variablecell = false, 
                   precond = nothing,
                   method = :auto,
                   g_tol = 1e-6, 
                   f_tol = 0.0, 
                   # should sub tolerances for units convert ... 
                  #  frc_tol = 1e-6u"eV/Å", 
                  #  vir_tol = 1e-6u"eV",
                  #  E_tol = 0.0u"eV", 
                   verbose = 1,
                   # robust_energy_difference = false, 
                   store_trace = false,
                   extended_trace = false,
                   maxstep = Inf,
                   callback = nothing,
                   g_calls_limit = 1_000, 
                   convert_sys = true )

   # convert the system to a format that supports setters 
   if convert_sys 
      sys1 = AosSystem(sys) 
   else 
      sys1 = sys 
   end 

   # create an objective function
   dofmgr = DofManager(sys1; variablecell = variablecell)
   x0 = get_dofs(sys, dofmgr)
   obj_f, obj_g! = get_obj_fg!(sys, calc, dofmgr)

   if isnothing(precond)
      precond = I 
   end 
   # # create a preconditioner
   # to be re-introduced 
   # if isa(precond, Symbol)
   #    if precond == :auto
   #       if fixedcell(at)
   #          precond = :exp
   #       else
   #          precond = :id
   #       end
   #    end
   #    if precond == :exp
   #       if method == :lbfgs
   #          precond = Exp(at, energyscale = :auto)
   #       else
   #          precond = Exp(at)
   #       end
   #    elseif precond == :id
   #       precond = I
   #    else
   #       error("unknown symbol for precond")
   #    end
   # end

   # choose the optimisation method Optim.jl
   if method == :auto || method == :cg
      if precond == I
         optimiser = ConjugateGradient(linesearch = BackTracking(order=2, maxstep=maxstep))
      else
         @assert precond isa AbstractMatrix
         optimiser = ConjugateGradient( P = precond,
                           precondprep = (P, x) -> precond, # update!(P, at, x),
                           linesearch = BackTracking(order=2, maxstep=maxstep) )
      end
   elseif method == :lbfgs
      optimiser = LBFGS( P = precond,
                        precondprep = (P, x) -> update!(P, at, x),
                        alphaguess = LineSearches.InitialHagerZhang(),
                        linesearch = BackTracking(order=2, maxstep=maxstep) )
   elseif method == :sd
      optimiser = Optim.GradientDescent( P = precond,
                  precondprep = (P, x) -> update!(P, at, x),
                  linesearch = BackTracking(order=2, maxstep=maxstep) )
   else
      error("GeomOpt.minimise : unknown `method` option")
   end

   results = optimize( obj_f, obj_g!, x0, optimiser,
                        Optim.Options( f_tol = f_tol, g_tol = g_tol,
                                       g_calls_limit = g_calls_limit,
                                       store_trace = store_trace,
                                       extended_trace = extended_trace,
                                       callback = callback,
                                       show_trace = (verbose > 1)) )
   set_dofs!(sys1, dofmgr, Optim.minimizer(results))
   # analyse the results
   if verbose > 0
      println(results)
   end

   return sys1, results
end
