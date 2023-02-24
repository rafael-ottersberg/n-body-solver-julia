# Install
Run julia, then:
```julia
using Pkg
Pkg.add(["CSV", "DataFrames", "BenchmarkTools"])
```

# Run
`julia Main.jl`

For interactive use:
Run julia, then:
```julia
using BenchmarkTools
include("NBodySolver.jl")
x, y, z, vx, vy, vz, m = readdata("solar_jfc.dat")
@benchmark leapfrog!(x, y, z, vx, vy, vz, m, 6.67408e-11, 24*60*60, 100)
```