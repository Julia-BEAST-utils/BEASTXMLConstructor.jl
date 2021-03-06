using Test, SafeTestsets


# @time @safetestset "Data simulation test" begin include("simulation_test.jl") end
# @time @safetestset "RTrees Test" begin include("RTrees_test.jl") end
@time @safetestset "XMLConstructor Test" begin include("xml_test.jl") end