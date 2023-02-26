using BenchmarkTools
include("NBodySolver.jl")
x, y, z, vx, vy, vz, m = readdata("solar_jfc.dat")
leapfrog!(x, y, z, vx, vy, vz, m, 6.67408e-11, 24*60*60, 100)
print(@benchmark leapfrog!(x, y, z, vx, vy, vz, m, 6.67408e-11, 24*60*60, 100))