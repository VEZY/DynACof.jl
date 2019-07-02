using DynACof
using Test

@test GDD(30.0,27.0,5.0,27.0)==0.0
@test GDD(25.,5.0,28.0)==20.0
@test GDD(5.0,5.0,28.0)==0.0

df= DataFrame(A = 1:10)
@test is_missing(df,"A")==false
@test is_missing(df,"B")==true

Parameters= Dict("Stocking_Coffee"=> 5580)
@test is_missing(Parameters,"Stocking_Coffee")==false
@test is_missing(Parameters,"B")==true
