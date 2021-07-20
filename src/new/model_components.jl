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

################################################################################
## TraitData
################################################################################


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

################################################################################
## Factor Model
################################################################################

struct FactorModel <: AbstractDataModel
    data::TraitData
    L::Matrix{Float64} # k X p
    prec::Vector{Float64} # factor precision
end

function FactorModel(data::TraitData, k::Int)
    p = size(data.data, 2)
    L = zeros(k, p)
    return FactorModel(data, L, ones(p))
end

function input_dim(model::FactorModel)
    return size(model.L, 1)
end

function model_elements(model::FactorModel;
            tree_model::GeneralizedXMLElement,
            is_submodel::Bool = false)
    trait_name = get_traitdata(model).trait_name
    id_header = is_submodel ? trait_name * "." : ""
    loadings = matrixParameterXML(model.L, id = id_header * "L")
    precision = parameterXML(id = id_header * "factorPrecision",
            value = model.prec,
            lower = zeros(length(model.prec)))


    ifa = integratedFactorModelXML(loadings = loadings,
        precision = precision,
        treeModel = tree_model,
        traitParameter = find_trait_parameter(tree_model, trait_name),
        trait_name = trait_name,
        id = id_header * ".factormodel")

    return Organizer([loadings, precision, ifa],
            likelihoods = [ifa],
            # priors = [], #TODO
            loggables = [loadings, precision],
            partial_providers = [ifa])

end


################################################################################
## MBD Model
################################################################################

struct RepeatedMeasuresModel <: AbstractDataModel
    data::TraitData
    # variance::Matrix{Float64}
end

function model_elements(model::RepeatedMeasuresModel;
            tree_model::GeneralizedXMLElement, is_submodel::Bool = false)
    # L_id = is_submodel ? get_traitdata(model).trait_name * ".L" : "L"
    # loadings_xml = matrixParameterXML(model.L, id=L_id)

    return Organizer(GeneralizedXMLElement[])

    # factor_model_id =
end


################################################################################
## Joint model
################################################################################


function decomposed_var_prior(var_mat::GeneralizedXMLElement, dim::Int)
    @assert var_mat.name == "compoundSymmetricMatrix"
    diag_param = find_element(var_mat, name = "diagonal", passthrough = true)
    offdiag_param = find_element(var_mat, name = "offDiagonal", passthrough = true)

    priors = [halfTPriorXML(diag_param), lkjPriorXML(offdiag_param, dim)]
    return Organizer(priors, priors = priors)
end

function make_xml(model::GeneralizedContinuousTraitModel)
    @unpack data, taxa, newick, models = model
    txxml = taxaXML(taxa, data)
    nxml = newickXML(newick)
    tmxml = treeModelXML(nxml, data, id = "treeModel")


    n_models = length(models)
    tree_dims = [input_dim(sub_model) for sub_model = models]
    q = sum(tree_dims)

    var_mat = compoundSymmetricMatrixXML(q, "variance.diagonal",
            "variance.offDiagonal", id="mbd.variance")
    # println(make_element(var_mat))
    p_mat = cachedMatrixInverseXML(var_mat, id = "mbd.precision")

    mbd_xml = mbdXML(p_mat, id="diffusionModel")

    elements = [txxml, nxml, tmxml, mbd_xml]
    org = Organizer(elements)

    for sub_model in models
        org = vcat(org,  model_elements(sub_model, tree_model = tmxml, is_submodel = true))
    end

    joint_extension = GeneralizedXMLElement("jointPartialsProvider",
            id = "jointModel",
            children = org.partial_providers)
    trait_likelihood = traitDataLikelihoodXML(
            diffusion_model = mbd_xml,
            tree_model = tmxml,
            extension_model = joint_extension,
            root_mean = zeros(q)
    )

    push!(org.elements, trait_likelihood)
    push!(org.likelihoods, trait_likelihood)

    var_priors = decomposed_var_prior(var_mat, q)
    org = vcat(org, var_priors)



    # TODO: likelihood
    # TODO: operators



    mcmc = mcmcXML(org, GeneralizedXMLElement("operators", id="TODO"),
            chain_length = 1000,
            screen_logEvery = 10,
            file_logEvery = 100,
            file_name = "test.log")



    beast = BEASTXMLDocument([org.elements; mcmc])

    return make_xml(beast)
end



