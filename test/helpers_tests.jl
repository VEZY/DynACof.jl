using DynACof
using Test
using DataFrames

@testset "GDD()" begin
    @test GDD(30.0,27.0,5.0,27.0)==0.0
    @test GDD(25.,5.0,28.0)==20.0
    @test GDD(5.0,5.0,28.0)==0.0
end;


@testset "warn_var()" begin
    @test_throws ErrorException warn_var("Date","Start_Date from Parameters","error")
    @test_throws ErrorException warn_var("Date")
end;


@testset "is_missing()" begin
    #For a dataframe:
    df= DataFrame(A = 1:10)
    @test is_missing(df,"A")==false
    @test is_missing(df,"B")==true

    #For a dictionary:
    Parameters= Dict("Stocking_Coffee"=> 5580)
    @test is_missing(Parameters,"Stocking_Coffee")==false
    @test is_missing(Parameters,"B")==true
end;

@testset "helpers" begin
    @test struct_to_tuple(constants, constants()) == (cp = 0.0010130000000000007, epsi = 0.622, pressure0 = 101.325, FPAR = 0.5, g = 9.81, Rd = 287.0586, Rgas = 8.314, Kelvin = 273.15, vonkarman = 0.41, MJ_to_W = 1.0000000000000006e-6, Gsc = 1367.0, σ = 5.670367e-8, H2OMW = 0.018, W_umol = 4.57, λ = 2.45, cl = 0.4, Dheat = 2.15e-5)
end;


@testset "ecophysiology helpers" begin
    @test rH_to_VPD(0.5,20.0) ≈ 1.1662980110489036
    @test esat(20.0,"Allen_1998") ≈ 2.3382812709274456
    @test esat_slope(20.0,"Allen_1998") ≈ 0.1447462277835135
    @test VPD_to_e(1.5, 25.0, "Sonntag_1990") ≈ 1.6600569164883336
    @test virtual_temp(20.0,1010.0,1.5) ≈ 20.091375550353973
    @test dew_point(20.0, 2.0) ≈ 11.252745464952273
end;

@testset "Meteorology helpers" begin
    @test cos°(0.0) == 1.0
    @test cos°(45.0) ≈ 0.707106781
    @test sin°(0.0) == 0.0
    @test sin°(90.0) == 1.0
    @test tan°(0.0) == 0.0
    @test tan°(45.0) == 1.0
    @test asin°(0.0) == 0.0
    @test asin°(1.0) == 90.0
    @test acos°(0.0) == 90.0
    @test acos°(1.0) == 0.0
    @test atan°(0.0) == 0.0
    @test atan°(1.0) == 45.0
    @test pressure_from_elevation(600.0, 25.0, 1.5) ≈ 94.63400648527185
    @test Rad_ext(1,35.0) ≈ 16.907304429496513
    @test Rad_ext(182,35.0) ≈ 41.482922899233984
    @test diffuse_fraction(1,25.0,35.0) == 0.23
    @test diffuse_fraction(1,1.0,35.0) == 1.0
    @test diffuse_fraction(1,10.0,35.0) ≈ 0.4664679
    @test sun_zenithal_angle(198,43.6109200) ≈ 0.3867544130665691 # considering solar time at noon
    @test Rad_net(1,5.0,16.0,10.0,1.5,9.0,1000.0,0.146) ≈ 4.2700530171252735
    @test days_without_rain([0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.2,0.6]) == [1,2,3,4,0,1,2,0,0]
end;