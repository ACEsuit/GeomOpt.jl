# GeomOpt

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ACEsuit.github.io/GeomOpt.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ACEsuit.github.io/GeomOpt.jl/dev/)
[![Build Status](https://github.com/ACEsuit/GeomOpt.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ACEsuit/GeomOpt.jl/actions/workflows/CI.yml?query=branch%3Amain)

**CLOSE BUT NOT QUITE READY FOR USE. PLEASE WAIT**

This is a straightforward translation of the geometry optimization interface in JuLIP.jl to the AtomsBase / AtomsCalculators eco-system. It is intended to be entirely temporary until [`GeometryOptimization.jl`](https://github.com/JuliaMolSim/GeometryOptimization.jl) matures sufficiently. 

At the moment of writing this package, it is in some ways more restrictive than GeometryOptimization and much more opinionated. But on the other hand it has a more structured approach to non-dimensionalization, it is easier to use, and aims to provide usable interfaces for saddle search and MEPs.

Because this is intended as a temporary package rather than an alternative to `GeometryOptimization.jl` I am currently not planning to register it in `General`, but only in the [ACE registry](https://github.com/ACEsuit/ACEregistry). To install the package first setup the ACE registry 
```
] registry add https://github.com/ACEsuit/ACEregistry.git
```
and then add the package as usual 
```
] add GeomOpt 
```

