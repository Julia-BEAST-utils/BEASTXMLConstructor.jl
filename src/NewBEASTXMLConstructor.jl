module NewBEASTXMLConstructor


export MBDModel,
       FactorModel,
       GeneralizedContinuousTraitModel,
       TraitData,
       make_xml,
       run_test


using LightXML, UnPack, LinearAlgebra


new_dir = "new"
include(joinpath(new_dir, "generalized.jl"))
include(joinpath(new_dir, "model_components.jl"))


end