using DynACof
using Test

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
end;
