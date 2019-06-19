using DynACof
using Test

@test GDD(30.0,27.0,5.0,27.0)==0.0
@test GDD(25.,5.0,28.0)==20.0
@test GDD(5.0,5.0,28.0)==0.0
