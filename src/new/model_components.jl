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


struct JointTraitModel{T <: AbstractDataModel, S <: AbstractMatrix{Float64}}
    data::DataPairs
    taxa::Vector{String}
    models::Vector{T}
    diffusion_variance::S
    newick::String
end

function JointTraitModel(models::Vector{<:AbstractDataModel},
        newick::String;
        V::AbstractMatrix{Float64} = Diagonal(ones(sum(input_dim.(models)))))
    n = length(models)

    #put factor models first
    priorities = priority.(models)
    perm = sortperm(priorities)
    models = models[perm]

    @assert n > 0
    taxa = get_taxa(models[1])
    for i = 2:n
        model_taxa = get_taxa(models[i])
        @assert model_taxa == taxa
    end

    data = [get_traitname(models[i]) => get_data(models[i]) for i = 1:n]
    return JointTraitModel(data, taxa, models, V, newick)
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

function TraitData(old_data::Matrix{Union{Missing, Float64}}, args...; kw_args...)
    data = zeros(size(old_data))
    for i = 1:length(data)
        if ismissing(old_data[i])
            data[i] = NaN
        else
            data[i] = old_data[i]
        end
    end
    return TraitData(data, args...; kw_args...)
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

function parse_traitdata(path::String; kwargs...)
    df = CSV.read(path, DataFrame)
    parse_traitdata(df; kwargs...)
end

function parse_traitdata(df::DataFrame; trait_name::String = "trait")
    nms = names(df)
    @assert nms[1] == "taxon"
    taxa = df.taxon
    data = Matrix(df[!, 2:end])
    trait_names = nms[2:end]

    return TraitData(data, taxa, trait_names, trait_name)
end


################################################################################
## Factor Model
################################################################################

struct FactorModel <: AbstractDataModel
    data::TraitData
    L::Matrix{Float64} # k X p
    prec::Vector{Float64} # factor precision
    standardize::Bool
end

function FactorModel(data::TraitData, k::Int; standardize::Bool = true)
    p = size(data.data, 2)
    L = zeros(k, p)
    return FactorModel(data, L, ones(p), standardize)
end

function input_dim(model::FactorModel)
    return size(model.L, 1)
end

function priority(::FactorModel)
    return 1
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
        id = id_header * ".factormodel",
        standardize = model.standardize)

    loadings_prior = distributionLikelihoodXML(data_parameter = loadings,
            distribution_model = standardNormalDistributionXML())
    precision_prior = gammaPriorXML(precision)


    return Organizer([loadings, precision, ifa, loadings_prior, precision_prior],
            likelihoods = [ifa],
            priors = [loadings_prior, precision_prior],
            loggables = [loadings, precision],
            partial_providers = [ifa])

end



function correlation_offset(offset::Int, dim::Int)
    corr_offset = 0
    for i = 1:offset #TODO: actually do math and replace for loop
        corr_offset += dim - i
    end
    return corr_offset
end

function set_correlation_mask_full!(x::AbstractVector{<:Real},
        model::FactorModel, dim::Int, offset::Int)

    k = input_dim(model)

    if offset == 0
        for row = 1:k
            for col = 1:dim
                x[(row - 1) * dim + col] = 0
            end
        end
    end


    for row = (offset + 1):(offset + k)
        for col = (offset + k + 1):dim
            x[dim * (row - 1) + col] = 0
        end
    end
end


function set_correlation_mask!(x::AbstractVector{<:Real}, model::FactorModel,
        dim::Int, offset::Int)
    k = input_dim(model)

    corr_offset = correlation_offset(offset, dim)

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
end

function setup_operators(::FactorModel, org::Organizer;
        trait_likelihood::GeneralizedXMLElement,
        loadings_operator::String = "hmc")

    @assert length(org.partial_providers) == 1
    provider = org.partial_providers[1]
    super_provider = find_partials_provider(trait_likelihood)

    loadings = find_loadings(provider)
    loadings_prior = find_prior(org, loadings)

    if lowercase(loadings_operator) == "hmc"

        super_model = super_provider !== provider

        loadings_gradient_providers = [provider, trait_likelihood]
        if super_model
            push!(loadings_gradient_providers, super_provider)
        end


        loadings_like_grad = GeneralizedXMLElement(
                "integratedFactorAnalysisLoadingsGradient",
                children = loadings_gradient_providers)



        loadings_prior_grad = gradientXML([loadings_prior, loadings])
        loadings_grad = jointGradientXML([loadings_like_grad, loadings_prior_grad])
        loadings_op = hmcXML(gradient = loadings_grad, parameter = loadings)
    elseif lowercase(loadings_operator) == "gibbs"
        loadings_op = GeneralizedXMLElement(bn.LOADINGS_GIBBS_OP,
                children = [provider, trait_likelihood, loadings_prior],
                attributes = [bn.WEIGHT => 1.0, bn.RANDOM_SCAN => false,
                              bn.NEW_MODE => true, bn.CONSTRAINT => "none",
                              bn.SPARSITY => "none"])
    else
        error("unrecognized loadings operator type '$loadings_operator'")
    end

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
    standardize::Bool

    function RepeatedMeasuresModel(data::TraitData, precision::AbstractMatrix{Float64};
                                   standardize::Bool = false)
        if !isposdef(precision)
            throw(ArgumentError("precision matrix must be positive-definite"))
        end

        return new(data, Matrix(precision), standardize)
    end
end

function RepeatedMeasuresModel(data::TraitData; kwargs...)
    p = size(data, 2)
    return RepeatedMeasuresModel(data, Matrix(Diagonal(ones(p))); kwargs...)
end

function priority(::RepeatedMeasuresModel)
    return 2
end

function model_elements(model::RepeatedMeasuresModel;
            tree_model::GeneralizedXMLElement, is_submodel::Bool = false)
    prec_id = is_submodel ?
            get_traitdata(model).trait_name * ".residualPrecision" :
                    "residualPrecision"
    prec = matrixParameterXML(model.precision, id = prec_id)
    trait_name = get_traitdata(model).trait_name

    rm = repeatedMeasuresModelXML(tree_model = tree_model, precision = prec,
            trait_name = trait_name,
            standardize = model.standardize)


    return Organizer([rm], loggables = [prec], partial_providers = [rm])

    # factor_model_id =
end

function setup_operators(::RepeatedMeasuresModel, org::Organizer;
        trait_likelihood::GeneralizedXMLElement, args...)
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


function set_correlation_mask_full!(x::AbstractVector{<:Real},
    model::RepeatedMeasuresModel, dim::Int, offset::Int)

    k = input_dim(model)
    for row = (offset + 1):(offset + k)
        for col = (row + 1):dim
            x[dim * (row - 1) + col] = 0
        end
    end

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

# function make_xml(model::JointTraitModel)

#     models = model.models
#     if length(models) > 1 && models[2] <: FactorModel
#         return make_multiFactor_xml(model)
#     else
#         return make_simple_joint_xml(model)
#     end
# end

function setup_taxa_and_tree(data::DataPairs, taxa::Vector{String}, newick::String)
    txxml = taxaXML(taxa, data)
    nxml = newickXML(newick)
    tmxml = treeModelXML(nxml, data, id = "treeModel")

    return (taxa_xml = txxml, newick_xml = nxml, treeModel_xml = tmxml)
end

import LinearAlgebra.cholesky
function cholesky(X::Hermitian{Float64, Diagonal{Float64, Vector{Float64}}})
    return cholesky(X.data)
end


function multiFactor_correlation_parameters(C::AbstractMatrix{Float64})
    n = size(C, 1)
    L = cholesky(C).L #need reverse cholesky (for now)
    decomp = matrixParameterXML(L, id="corr.decomp")
    inner_prod = transformedParameterXML(decomp,
        PassthroughXMLElement("matrixInnerProductTransform", decomp),
        as_matrix=true,
        id="corr.full")
    mask = ones(Int, n, n)
    for i = 1:n
        for j = 1:(i - 1)
            mask[i, j] = 0
        end
    end

    masked = maskedParameterXML(inner_prod, vec(mask'), id="corr.upper")

    return (decomposed_corr = decomp, full_corr = inner_prod, masked_corr = masked)
end

function multiFactor_correlation_parameters(n::Int)
    multiFactor_correlation_parameters(Diagonal(ones(n)))
end




function make_xml(model::JointTraitModel;
        mcmc_options = MCMCOptions(),
        file_name::String = "test.log",
        loadings_operator::String = "hmc")
    @unpack data, taxa, newick, models = model

    is_multiFactor = count(x -> typeof(x) <: FactorModel, models) > 1


    taxa_and_tree = setup_taxa_and_tree(data, taxa, newick)
    org = Organizer(collect(taxa_and_tree))
    @unpack taxa_xml, newick_xml, treeModel_xml = taxa_and_tree

    n_models = length(models)
    tree_dims = [input_dim(sub_model) for sub_model = models]
    q = sum(tree_dims)


    if is_multiFactor
        S = inv(Diagonal(sqrt.(diag(model.diffusion_variance))))
        C = Hermitian(S * model.diffusion_variance * S)

        corr_params = multiFactor_correlation_parameters(C)
        @unpack decomposed_corr, full_corr, masked_corr = corr_params
        decomposed_corr_prior = determinantPriorXML(decomposed_corr)


        push!(org, decomposed_corr, loggable = true)
        push!(org, decomposed_corr_prior, prior = true)
        # push!(org, full_corr, loggable=true)


        diag_var = parameterXML(value=diag(model.diffusion_variance), lower = zeros(q), id="variance.diagonal")

        var_mat = compoundSymmetricMatrixXML(diag_var, masked_corr,
                as_correlation = true, is_cholesky = false,
                strictly_upper = false, id = "variance")
    else
        var_mat = compoundSymmetricMatrixXML(q, "variance.diagonal",
                "variance.offDiagonal", id="mbd.variance")
    end
    p_mat = cachedMatrixInverseXML(var_mat, id = "mbd.precision")

    mbd_xml = mbdXML(p_mat, id="diffusionModel")

    push!(org, var_mat, loggable=true)
    push!(org, mbd_xml)

    sub_model_components = Vector{Organizer}(undef, n_models)

    for i = 1:n_models
        sub_model = models[i]
        components = model_elements(sub_model, tree_model = treeModel_xml, is_submodel = true)
        sub_model_components[i] = components
        org = vcat(org,  components)
    end

    joint_extension = GeneralizedXMLElement("jointPartialsProvider",
            id = "jointModel",
            children = org.partial_providers)
    trait_likelihood = traitDataLikelihoodXML(
            diffusion_model = mbd_xml,
            tree_model = treeModel_xml,
            extension_model = joint_extension,
            root_mean = zeros(q)
    )

    push!(org, trait_likelihood, likelihood = true)

    # diffusion priors + operator
    operators = GeneralizedXMLElement[]


    if is_multiFactor
        diff_like_grad = diffusionGradientXML(trait_likelihood = trait_likelihood,
                precision_parameter = p_mat,
                id = "corr.likelihood.gradient",
                parameter = "decomposedCorrelation")

        orthogonality_structure = Vector{Int}[]
        ind = 1
        for i = 1:n_models
            next_ind = ind + tree_dims[i]
            if typeof(models[i]) <: FactorModel
                rows = collect(ind:(next_ind - 1))
                Base.push!(orthogonality_structure, rows)
            end
            ind = next_ind
        end

        corr_mask = ones(Int, q * q)
        offset = 0
        for sub_model in models
            set_correlation_mask_full!(corr_mask, sub_model, q, offset)
            offset += input_dim(sub_model)
        end

        prior_grad = GeneralizedXMLElement("gradient",
            children = [decomposed_corr_prior, decomposed_corr])

        hmc_op = hmcXML(
                gradient = GeneralizedXMLElement("jointGradient",
                        children = [diff_like_grad, prior_grad]),
                parameter = decomposed_corr,
                is_geodesic = true,
                orthogonality_structure = orthogonality_structure
                # mask_parameter = parameterXML(value = corr_mask)
                )
        push!(operators, hmc_op)
    else
        diag_param = find_element(var_mat, name = "diagonal", passthrough = true)
        offdiag_param = find_element(var_mat, name = "offDiagonal", passthrough = true)

        diag_prior = halfTPriorXML(diag_param)
        offdiag_prior = lkjPriorXML(offdiag_param, q)

        var_priors = [diag_prior, offdiag_prior]
        push!.(Ref(org), var_priors, prior = true)

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

        push!(operators, hmc_op)
    end

    # rw_mask = diff_mask[(q + 1):end] #should create new array w/ unlinked elements
    # for i = 1:length(rw_mask)
    #     rw_mask[i] = rw_mask[i] == 0 ? 1 : 0
    # end
    # rw_op = randomWalkXML(maskedParameterXML(offdiag_param, corr_mask),
    #     window_size = 0.1)

    for i = 1:n_models
        sub_model = models[i]
        sub_components = sub_model_components[i]
        sub_operators = setup_operators(sub_model, sub_components,
                trait_likelihood = trait_likelihood,
                loadings_operator = loadings_operator)
        operators = [operators; sub_operators]
    end


    ops = GeneralizedXMLElement("operators", children = operators, id = "operators")



    push!(org, ops)
    push!(org, GeneralizedXMLElement("correlationMatrix",
            id = "correlation",
            child = var_mat),
            loggable=true)
    push!(org, traitLoggerXML(tree_model = treeModel_xml,
            trait_likelihood = trait_likelihood),
            loggable=true)

    mcmc = mcmcXML(org, ops, mcmc_options, file_name = file_name)



    beast = BEASTXMLDocument([org.elements; mcmc])

    return make_xml(beast)
end

function save_xml(model::JointTraitModel, path::String; kw_args...)
    nm = basename(path)
    if endswith(path, ".xml")
        nm = nm[1:(end - 3)] * "log"
    end

    xml = make_xml(model, file_name = nm; kw_args...)

    LightXML.save_file(xml, path)
    free(xml)
    return nothing
end



