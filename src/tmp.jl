using Revise
using DynACof

file= "D:/OneDrive/Travail_These/Test_DCM/1-Input/Aquiares/Meteorology.txt"
using CSV
MetData= CSV.read(file; copycols=true)