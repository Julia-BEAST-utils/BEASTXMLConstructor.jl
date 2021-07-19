module NewBEASTXMLConstructor


export MBDModel,
       FactorModel,
       GeneralizedContinuousTraitModel,
       TraitData,
       make_xml


using LightXML, UnPack, LinearAlgebra


new_dir = "new"
include(joinpath(new_dir, "generalized.jl"))

end