using Pkg
Pkg.activate(".")
Pkg.instantiate()

using Revise

includet("q1.jl")
includet("q2.jl")
includet("q3.jl")
includet("q3_tfunctions.jl")
includet("q3_tests.jl")

println("New setup complete")