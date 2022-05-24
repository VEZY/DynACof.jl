file = download("https://raw.githubusercontent.com/VEZY/DynACof.jl_inputs/master/meteorology.txt")

@testset "Meteo computations" begin
    p = import_parameters("package")
    m = meteorology(file, p)

    @test size(m) == (13880, 19)
    @test [m.Rain[1], m.Rain[end]] == [5.29, 0.0]
    @test m.DaysWithoutRain[1:2] == [0.0, 1.0]
    @test maximum(m.DaysWithoutRain) == 62.0

    # Using sums as they integrate the results for the whole sequence:
    @test sum(m.DaysWithoutRain) == 17152.0
    @test sum(m.DegreeDays) ≈ 133421.951 atol = 1e-2
    @test sum(m.FDiff) ≈ 5300.898 atol = 1e-2
    @test sum(m.Rn) ≈ 179314.01177 atol = 1e-2
    @test sum(m.ZEN) ≈ 3880.577 atol = 1e-2
end;


@testset "Meteo with sub-period" begin
    p = import_parameters("package")
    m = meteorology(file, p, ["2015-01-01", "2015-12-31"])

    @test size(m) == (365, 19)
    @test m.Date[1] == Dates.Date("2015-01-01")
    @test m.Date[end] == Dates.Date("2015-12-31")
end
