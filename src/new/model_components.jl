
################################################################################
## Model Structures
################################################################################

abstract type AbstractDataModel end

function output_dim(model::AbstractDataModel)
    return size(get_traitdata(model), 1)
end

function input_dim(model::AbstractDataModel)
    return output_dim(model)
end


struct GeneralizedContinuousTraitModel{T <: AbstractDataModel}
    data::DataPairs
    taxa::Vector{String}
    models::Vector{T}
    newick::String

end

function GeneralizedContinuousTraitModel(models::Vector{<:AbstractDataModel},
        newick::String)
    n = length(models)
    @assert n > 0
    taxa = get_taxa(models[1])
    for i = 2:n
        model_taxa = get_taxa(models[i])
        @assert model_taxa == taxa
    end

    data = [get_traitname(models[i]) => get_data(models[i]) for i = 1:n]
    return GeneralizedContinuousTraitModel(data, taxa, models, newick)
end


struct TraitData
    data::Matrix{Float64}
    taxa::Vector{String}
    trait_names::Vector{String}
    trait_name::String
end

import Base.size
function size(td::TraitData, args...)
    return size(td.data, args...)
end

function TraitData(data::Matrix{Float64}, taxa::Vector{String};
        trait_name::String = "traits",
        trait_names::Vector{String} = ["trait_$i" for i = 1:size(data, 2)])
    return TraitData(data, taxa, trait_names, trait_name)
end

function get_taxa(model::AbstractDataModel)
    return get_traitdata(model).taxa
end

function get_traitname(model::AbstractDataModel)
    return get_traitdata(model).trait_name
end

function get_data(model::AbstractDataModel)
    return get_traitdata(model).data
end


function get_traitdata(model::AbstractDataModel)::TraitData
    return model.data
end



struct FactorModel <: AbstractDataModel
    data::TraitData
    L::Matrix{Float64} # k X p
end

function FactorModel(data::TraitData, k::Int)
    p = size(data.data, 2)
    L = zeros(k, p)
    return FactorModel(data, L)
end

function input_dim(model::FactorModel)
    return size(model.L, 1)
end

function model_elements(model::FactorModel, ; is_submodel::Bool = false)
    L_id = is_submodel ? get_traitdata(model).trait_name * ".L" : "L"
    loadings_xml = matrixParameterXML(model.L, id=L_id)

    # factor_model_id =
end


struct MBDModel <: AbstractDataModel
    data::TraitData
    # variance::Matrix{Float64}
end



function make_xml(model::GeneralizedContinuousTraitModel)
    @unpack data, taxa, newick, models = model
    txxml = taxaXML(taxa, data)
    nxml = newickXML(newick)
    tmxml = treeModelXML(nxml, data)

    n_models = length(models)
    tree_dims = [input_dim(sub_model) for sub_model = models]
    q = sum(tree_dims)

    var_mat = compoundSymmetricMatrixXML(q, "variance.diagonal",
            "variance.offDiagonal", id="mbd.variance")
    # println(make_element(var_mat))
    p_mat = cachedMatrixInverseXML(var_mat, id = "mbd.precision")

    mbd_xml = mbdXML(p_mat)


    beast = BEASTXMLDocument([txxml, nxml, tmxml, mbd_xml])

    return make_xml(beast)
end



