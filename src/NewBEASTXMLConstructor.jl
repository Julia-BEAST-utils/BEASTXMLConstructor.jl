module NewBEASTXMLConstructor


export MBDModel,
       FactorModel,
       GeneralizedContinuousTraitModel,
       TraitData,
       make_xml,
       run_test,
       FactorModel,
       RepeatedMeasuresModel,
       AbstractDataModel


using LightXML, UnPack, LinearAlgebra

include("BeastNames.jl")
const bn = BeastNames

const DataPairs = Vector{Pair{String, Matrix{Float64}}}
const Nullable{T} = Union{Nothing, T} where T <: Any

new_dir = "new"
include(joinpath(new_dir, "generalized.jl"))
include(joinpath(new_dir, "wrappers.jl"))
include(joinpath(new_dir, "model_components.jl"))
include(joinpath(new_dir, "elements.jl"))




function run_test()
    n = 10
    taxa = ["taxon_$i" for i = 1:n]


    ps = [3, 4, 5]
    nms = ["traitA", "traitB", "traitC"]
    data = [TraitData(randn(n, ps[i]), taxa, trait_name = nms[i]) for i = 1:length(ps)]

    fac1 = FactorModel(data[1], 1)
    rm = RepeatedMeasuresModel(data[2])
    fac2 = FactorModel(data[3], 3)

    newick = "(((taxon_1:0.14391733598401604,(taxon_2:1.8282365833429757,taxon_8:0.2415088669042653):1.6685177420752977):1.9757457480210345,((taxon_4:0.12335520066265246,((taxon_7:1.037431406660496,taxon_6:1.188747348685923):0.2541982687655974,taxon_10:1.7618786450343555):0.8564170824143821):0.17352743115993857,(taxon_9:0.533959250714578,taxon_3:1.9969429483252148):0.0908838854524467):0.06602845267356403):0.6112158700724324,taxon_5:0.10024590507839154);"

    model = GeneralizedContinuousTraitModel([fac1, rm, fac2], newick)

    xml = make_xml(model)
    save_file(xml, raw"C:\Users\gabeh\OneDrive\Desktop\test.xml")
    print(xml)
    free(xml)
    return nothing
end

end