using DynACof
using Test

@testset "Test on helpers" begin
   include("helpers_tests.jl")
end

@testset "Tests on conductances" begin
   include("conductances_tests.jl")
end

@testset "Tests on diseases" begin
   include("disease_tests.jl")
end


@testset "parameters" begin
    include("parameters.jl")
end
