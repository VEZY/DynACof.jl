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


@testset "ecophysiology helpers" begin
    @test rH_to_VPD(0.5,20.0) ≈ 1.1662980110489036
    @test esat(20.0,"Allen_1998") ≈ 2.3382812709274456
    @test esat_slope(20.0,"Allen_1998") ≈ 0.1447462277835135
    @test VPD_to_e(1.5, 25.0, formula= "Sonntag_1990") ≈ 1.6600569164883336
    @test virtual_temp(20.0,1010.0,1.5) ≈ 20.091375550353973
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
end;
