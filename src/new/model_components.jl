################################################################################
## Model Structures
################################################################################

abstract type AbstractDataModel end

function output_dim(model::AbstractDataModel)
    return size(get_traitdata(model), 2)
end

function input_dim(model::AbstractDataModel)
    return output_dim(model)
end


struct JointTraitModel{T <: AbstractDataModel}
    data::DataPairs
    taxa::Vector{String}
    models::Vector{T}
    newick::String
end

function JointTraitModel(models::Vector{<:AbstractDataModel},
        newick::String)
    n = length(models)
    @assert n > 0
    taxa = get_taxa(models[1])
    for i = 2:n
        model_taxa = get_taxa(models[i])
        @assert model_taxa == taxa
    end

    data = [get_traitname(models[i]) => get_data(models[i]) for i = 1:n]
    return JointTraitModel(data, taxa, models, newick)
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

    loadings_prior = distributionLikelihoodXML(data_parameter = loadings,
            distribution_model = standardNormalDistributionXML())
    precision_prior = gammaPriorXML(precision)


    return Organizer([loadings, precision, ifa, loadings_prior, precision_prior],
            likelihoods = [ifa],
            priors = [loadings_prior, precision_prior],
            loggables = [loadings, precision],
            partial_providers = [ifa])

end

function set_correlation_mask!(x::AbstractVector{<:Real}, model::FactorModel,
        dim::Int, offset::Int)
    k = input_dim(model)

    corr_offset = 0
    for i = 1:offset #TODO: actually do math and replace for loop
        corr_offset += dim - i
    end

    for row = (offset + 1):(k + offset - 1)
        for col = 1:(offset - row + k)
            x[corr_offset + col] = 0
        end
    end
    return nothing
end



function set_diffusion_mask!(x::Vector{<:Real}, model::FactorModel,
        dim::Int, offset::Int)
    k = input_dim(model)
    big_k = offset + k
    x[(offset + 1):big_k] .= 0

    off_diagonals = @view x[(dim + 1):end]

    corr_offset = 0
    for row = 1:(big_k - 1)
        for col = max(1, offset - row + 1):(offset - row + k)
            off_diagonals[corr_offset + col] = 0
        end
        corr_offset += dim - row
    end

    # TODO: remove below
    @show offset

    L = zeros(dim, dim)
    ind = 0
    for i = 1:dim
        L[i, i] = x[i]
        for j = (i + 1):dim
            ind += 1
            L[i, j] = off_diagonals[ind]
            L[j, i] = L[i, j]
        end
    end
    display(L)
end

function setup_operators(::FactorModel, org::Organizer;
        trait_likelihood::GeneralizedXMLElement)

    @assert length(org.partial_providers) == 1
    provider = org.partial_providers[1]
    super_provider = find_partials_provider(trait_likelihood)

    super_model = super_provider !== provider

    loadings_gradient_providers = [provider, trait_likelihood]
    if super_model
        push!(loadings_gradient_providers, super_provider)
    end

    loadings_like_grad = GeneralizedXMLElement(
            "integratedFactorAnalysisLoadingsGradient",
            children = loadings_gradient_providers)

    loadings = find_loadings(provider)
    loadings_prior = find_prior(org, loadings)

    loadings_prior_grad = gradientXML([loadings_prior, loadings])
    loadings_grad = jointGradientXML([loadings_like_grad, loadings_prior_grad])
    loadings_op = hmcXML(gradient = loadings_grad, parameter = loadings)

    precision = find_factor_precision(provider)
    prec_provider = normalExtensionXML(extension = provider, likelihood = trait_likelihood)
    prec_prior = find_prior(org, precision)
    prec_op = normalGammaGibbsXML(prior = prec_prior, provider = prec_provider)
    return [loadings_op, prec_op]

end

function find_prior(org::Organizer, parameter::GeneralizedXMLElement)
    id = prior_id(parameter)

    for prior in org.priors
        if prior.id == id
            return prior
        end
    end
    error("could not locate prior")
end



################################################################################
## MBD Model
################################################################################

struct RepeatedMeasuresModel <: AbstractDataModel
    data::TraitData
    precision::Matrix{Float64}
end

function RepeatedMeasuresModel(data::TraitData)
    p = size(data, 2)
    return RepeatedMeasuresModel(data, Matrix(Diagonal(ones(p))))
end

function model_elements(model::RepeatedMeasuresModel;
            tree_model::GeneralizedXMLElement, is_submodel::Bool = false)
    prec_id = is_submodel ?
            get_traitdata(model).trait_name * ".residualPrecision" :
                    "residualPrecision"
    prec = matrixParameterXML(model.precision, id = prec_id)
    trait_name = get_traitdata(model).trait_name

    rm = repeatedMeasuresModelXML(tree_model = tree_model, precision = prec,
            trait_name = trait_name)


    return Organizer([rm], loggables = [prec], partial_providers = [rm])

    # factor_model_id =
end

function setup_operators(::RepeatedMeasuresModel, org::Organizer;
        trait_likelihood::GeneralizedXMLElement)
    #TODO
    return GeneralizedXMLElement[]
end
function set_diffusion_mask!(::Vector{<:Real}, ::RepeatedMeasuresModel,
        ::Int, ::Int)
    # do nothing
end

function set_correlation_mask!(::Vector{<:Real}, ::RepeatedMeasuresModel,
    ::Int, ::Int)
# do nothing
end


################################################################################
## Joint model
################################################################################


# function decomposed_var_prior(var_mat::GeneralizedXMLElement, dim::Int)
#     @assert var_mat.name == "compoundSymmetricMatrix"
#     diag_param = find_element(var_mat, name = "diagonal", passthrough = true)
#     offdiag_param = find_element(var_mat, name = "offDiagonal", passthrough = true)

#     priors = [halfTPriorXML(diag_param), lkjPriorXML(offdiag_param, dim)]
#     return Organizer(priors, priors = priors)
# end

# function variance_gradients(;precision_parameter::GeneralizedXMLElement,
#             trait_lieklhiood::GeneralizedXMLElement,
#             diagonal_prior::GeneralizedXMLElement,
#             offdiagonal_prior::GeneralizedXMLElement,
#             )

function make_xml(model::JointTraitModel)
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
    org = Organizer(elements, loggables = [var_mat])

    sub_model_components = Vector{Organizer}(undef, n_models)

    for i = 1:n_models
        sub_model = models[i]
        components = model_elements(sub_model, tree_model = tmxml, is_submodel = true)
        sub_model_components[i] = components
        org = vcat(org,  components)
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

    # diffusion priors + operator

    diag_param = find_element(var_mat, name = "diagonal", passthrough = true)
    offdiag_param = find_element(var_mat, name = "offDiagonal", passthrough = true)

    diag_prior = halfTPriorXML(diag_param)
    offdiag_prior = lkjPriorXML(offdiag_param, q)

    var_priors = [diag_prior, offdiag_prior]


    diff_prior_grad = compoundGradientXML(
            [[diag_prior, diag_param], [offdiag_prior]],
            id="variance.prior.gradient")
    diff_like_grad = diffusionGradientXML(trait_likelihood = trait_likelihood,
            precision_parameter = p_mat,
            id = "variance.likelihood.gradient")

    diff_grad = jointGradientXML([diff_prior_grad, diff_like_grad],
            id = "variance.gradient")

    n_diff = div(q * (q + 1), 2)
    diff_mask = ones(Int, n_diff)
    corr_mask = ones(Int, n_diff - q)
    offset = 0
    for sub_model in models
        set_diffusion_mask!(diff_mask, sub_model, q, offset)
        set_correlation_mask!(corr_mask, sub_model, q, offset)
        offset += input_dim(sub_model)
    end


    diff_sub = @view diff_mask[(q + 1):end] #temporary solution
    for i = 1:length(corr_mask)
        if diff_sub[i] == 1 && corr_mask[i] == 1
            corr_mask[i] = 0
        end
    end

    diff_transform = multivariateCompoundTransformXML(
                transformXML("log", dim = q),
                lkjTransformXML(q)
    )


    hmc_op = hmcXML(gradient = diff_grad,
            parameter = GeneralizedXMLElement("compoundParameter",
                    children = [diag_param, offdiag_param]),
            mask_parameter = parameterXML(value = diff_mask),
            transform = diff_transform)

    # rw_mask = diff_mask[(q + 1):end] #should create new array w/ unlinked elements
    # for i = 1:length(rw_mask)
    #     rw_mask[i] = rw_mask[i] == 0 ? 1 : 0
    # end
    rw_op = randomWalkXML(maskedParameterXML(offdiag_param, corr_mask),
        window_size = 0.1)



    operators = GeneralizedXMLElement[hmc_op, rw_op]
    for i = 1:n_models
        sub_model = models[i]
        sub_components = sub_model_components[i]
        sub_operators = setup_operators(sub_model, sub_components,
                trait_likelihood = trait_likelihood)
        operators = [operators; sub_operators]
    end


    ops = GeneralizedXMLElement("operators", children = operators, id = "operators")




    org = vcat(org, Organizer(var_priors, priors = var_priors))
    push!(org.elements, ops)
    push!(org.loggables, GeneralizedXMLElement("correlationMatrix",
            id = "correlation",
            child = var_mat))
    push!(org.loggables, traitLoggerXML(tree_model = tmxml,
            trait_likelihood = trait_likelihood))





    # TODO: likelihood
    # TODO: operators



    mcmc = mcmcXML(org, ops,
            chain_length = 1000,
            screen_logEvery = 10,
            file_logEvery = 100,
            file_name = "test.log")



    beast = BEASTXMLDocument([org.elements; mcmc])

    return make_xml(beast)
end

function save_xml(model::JointTraitModel, path::String)
    xml = make_xml(model)
    LightXML.save_file(xml, path)
    free(xml)
end



