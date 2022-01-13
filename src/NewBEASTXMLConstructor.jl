module NewBEASTXMLConstructor


export MBDModel,
       FactorModel,
       JointTraitModel,
       TraitData,
       make_xml,
       run_test,
       FactorModel,
       RepeatedMeasuresModel,
       AbstractDataModel,
       save_xml,
       MCMCOptions


using LightXML, UnPack, LinearAlgebra

include("MCMCOptionsProvider.jl")
MCMCOptions = MCMCOptionsProvider.MCMCOptions

include("BeastNames.jl")
const bn = BeastNames

const DataPairs = Vector{Pair{String, Matrix{Float64}}}
const Nullable{T} = Union{Nothing, T} where T <: Any

new_dir = "new"
include(joinpath(new_dir, "generalized.jl"))
include(joinpath(new_dir, "wrappers.jl"))
include(joinpath(new_dir, "model_components.jl"))
include(joinpath(new_dir, "elements.jl"))



end