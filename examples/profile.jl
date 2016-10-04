using KSVD

include("benchmark.jl")


Y = rand(128, 100)

benchmark(ksvd, Y, 512)
