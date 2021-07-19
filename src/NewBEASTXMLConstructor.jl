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
include(joinpath(new_dir, "elements.jl"))




function run_test()
    n = 10
    taxa = ["taxon_$i" for i = 1:n]


    ps = [3, 5]
    nms = ["traitA", "traitB"]
    data = [TraitData(randn(n, ps[i]), taxa, trait_name = nms[i]) for i = 1:length(ps)]

    mbd = MBDModel(data[1])
    fac = FactorModel(data[2], 2)
    newick = "PUT NEWICK HERE"

    model = GeneralizedContinuousTraitModel([mbd, fac], newick)

    xml = make_xml(model)
    println(xml)
end

end