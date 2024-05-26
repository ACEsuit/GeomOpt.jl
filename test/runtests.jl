using GeomOpt
using Test

@testset "GeomOpt.jl" begin
    @testset "DofManager" begin include("test_dofmgr.jl") end 
    @testset "Optim" begin include("test_optim.jl") end 
end
