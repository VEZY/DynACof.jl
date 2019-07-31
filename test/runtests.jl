using DynACof
using Test

@testset "Test on helpers" begin
   include("helpers_tests.jl")
end

@testset "Tests on conductances" begin
   include("conductances_tests.jl")
end