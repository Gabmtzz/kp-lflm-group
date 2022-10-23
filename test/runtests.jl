using Test
using KPpack

@testset verbose = true "KPpack tests" begin
    include("testParams.jl")
    include("testSymb.jl")
end