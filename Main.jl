using BenchmarkTools
include("NBodySolver.jl")
bodies = readdata("solar_jfc.dat")
leapfrog!(bodies, 6.67408e-11, 24*60*60, 100)
trial = @benchmark leapfrog!(bodies, 6.67408e-11, 24*60*60, 100)
print(mean(trial))