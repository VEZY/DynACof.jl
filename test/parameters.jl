using DynACof
using Test

@testset "parameters" begin
    @test read_param_file(:constants,"package") == constants()
    a= import_parameters("package")
    @test length(a) == 208
    @test a.cp == constants().cp
end;
