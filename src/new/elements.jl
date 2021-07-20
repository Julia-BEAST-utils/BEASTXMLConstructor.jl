
################################################################################
## FactorModel
################################################################################

function integratedFactorModelXML(;
        loadings::GeneralizedXMLElement,
        precision::GeneralizedXMLElement,
        treeModel::GeneralizedXMLElement,
        traitParameter::GeneralizedXMLElement,
        trait_name::String,
        id::Union{String, Nothing} = nothing,
        standardize::Bool = true
        )
    children = [PassthroughXMLElement("loadings", loadings),
                PassthroughXMLElement("precision", precision),
                treeModel,
                PassthroughXMLElement("traitParameter", traitParameter)]
    attributes = [bn.TRAIT_NAME => trait_name, bn.STANDARDIZE => "$standardize"]
    return GeneralizedXMLElement("integratedFactorModel", id, attributes, children, nothing)
end

################################################################################
## Data
################################################################################

function taxaXML(taxa::Vector{String},
        data::DataPairs;
        dates::Union{Nothing, Vector{Float64}} = nothing)

    n = length(taxa)
    p = length(data)

    children = Vector{GeneralizedXMLElement}(undef, n)
    for i = 1:n
        date = isnothing(dates) ? nothing : dates[i]
        traits = [data[j][1] => data[j][2][i, :] for j = 1:p]
        children[i] = taxonXML(taxa[i], traits, date)
    end
    return GeneralizedXMLElement("taxa", children = children)
end

function taxonXML(taxon::String,
        traits::Vector{<:Pair{String, <:AbstractVector{Float64}}},
        date::Union{Float64, Nothing})
    children = [traitAttr(trait) for trait in traits]
    if !isnothing(date)
        error("not yet implemented")
    end
    return GeneralizedXMLElement("taxon", id = taxon, children = children)
end

function traitAttr(attr::Pair{String, <:AbstractArray{Float64}})
    return GeneralizedXMLElement("attr", attributes = ["name" => attr[1]],
        content = attr[2])
end

################################################################################
## Newick
################################################################################
function newickXML(newick::AbstractString;
        using_heights::Bool = true,
        id::String = "startingTree")
    using_dates = !using_heights

    attributes = ["usingHeights" => using_heights, "usingDates" => using_dates]
    return GeneralizedXMLElement("newick", attributes = attributes,
            content = newick, id = id)
end

################################################################################
## treeModel
################################################################################

function treeModelXML(tree::GeneralizedXMLElement, data::DataPairs;
        id::StringOrNothing = nothing)
    children = [tree,
            PassthroughXMLElement("rootHeight",
                parameterXML(id="treeModel.rootHeight")),
            nodeHeightsXML("treeModel.internalNodeHeights",
                internal_nodes = true, root_node = false),
            nodeHeightsXML("treeModel.allInternalNodeHeights",
                internal_nodes = true, root_node = true)]
    for i = 1:length(data)
        trait_name = data[i][1]
        trait_dim = size(data[i][2], 2)
        trait_xml = nodeTraitsXML(trait_name, trait_dim, "$trait_name.leafTraits")

        push!(children, trait_xml)
    end

    return GeneralizedXMLElement("treeModel", children = children, id = id)
end

function nodeHeightsXML(parameter::AbstractGeneralizedXMLElement;
        internal_nodes::Bool = true, root_node::Bool = true)
    attrs = ["internalNodes" => internal_nodes, "rootNode" => root_node]
    return GeneralizedXMLElement("nodeHeights",
        attributes = attrs,
        children=[parameter])
end

function nodeHeightsXML(parameter_id::String; kw_args...)
    return nodeHeightsXML(parameterXML(id=parameter_id); kw_args...)
end

function nodeTraitsXML(trait_name::String, dim::Int,
        parameter::AbstractGeneralizedXMLElement;
        internal_nodes::Bool = false, root_node::Bool = false,
        leaf_nodes::Bool = true,
        as_matrix::Bool = true)
    attrs = ["name" => trait_name,
             "traitDimension" => dim,
             "internalNodes" => internal_nodes,
             "rootNode" => root_node,
             "leafNodes" => leaf_nodes,
             "asMatrix" => as_matrix]
    return GeneralizedXMLElement("nodeTraits",
        attributes = attrs,
        children=[parameter])
end

function nodeTraitsXML(trait_name::String, dim::Int,
        parameter_id::AbstractString; kw_args...)
    return nodeTraitsXML(trait_name, dim, parameterXML(id=parameter_id);
            kw_args...)
end

function find_trait_parameter(tm::GeneralizedXMLElement, trait::String)
    @assert tm.name == "treeModel"
    node_traits = find_element(tm, attributes = Dict("name" => trait))

    return find_element(node_traits, name = "parameter")
end


################################################################################
## multivariateDiffusionModel
################################################################################

function mbdXML(precision_parameter::AbstractGeneralizedXMLElement;
        id::String = "diffusionModel")
    return GeneralizedXMLElement("multivariateDiffusionModel",
            children=[precision_parameter])
end


################################################################################
## parameters
################################################################################

const ArrayOrNothing{T} = Union{AbstractArray{T}, Nothing} where T
const PARAMETER_NAMES = ("value", "upper", "lower", "dimension")


function parameterXML(;id::String = nothing,
        value::ArrayOrNothing{Float64} = nothing,
        upper::ArrayOrNothing{Float64} = nothing,
        lower::ArrayOrNothing{Float64} = nothing,
        dimension::Union{Int, Nothing} = nothing)

    attributes = Pair{String, Any}[]
    values = (value, upper, lower, dimension)
    for i = 1:length(values)
        add_parameter_attribute!(attributes, PARAMETER_NAMES[i], values[i])
    end

    return GeneralizedXMLElement("parameter", id=id, attributes = attributes)
end

function add_parameter_attribute!(::Vector{<:Pair{String, <:Any}},
        ::String, ::Nothing)
    # do nothing
end

function add_parameter_attribute!(attributes::Vector{<:Pair{String, <:Any}},
        attr_name::String, attr_value::AbstractArray{Float64})
    push!(attributes, attr_name => attr_value)
end

function cachedMatrixInverseXML(parameter::AbstractGeneralizedXMLElement;
        id::String = nothing)
    return GeneralizedXMLElement("cachedMatrixInverse", id=id, child=parameter)
end

function compoundSymmetricMatrixXML(
        diagonal_parameter::AbstractGeneralizedXMLElement,
        offdiagonal_parameter::AbstractGeneralizedXMLElement;
        as_correlation::Bool = true,
        is_cholesky::Bool = true,
        id::StringOrNothing = nothing)

    attrs = ["asCorrelation" => as_correlation,
             "isCholesky" => is_cholesky]
    children = [PassthroughXMLElement("diagonal", diagonal_parameter),
                PassthroughXMLElement("offDiagonal", offdiagonal_parameter)]
    return GeneralizedXMLElement("compoundSymmetricMatrix", attributes = attrs,
            children = children, id = id)
end

function compoundSymmetricMatrixXML(X::AbstractArray{Float64, 2},
        diag_id::String, offdiag_id::String; is_posdef::Bool = true, kw_args...)
    @assert issymmetric(X)
    n = size(X, 1)

    lower_bounds = is_posdef ? zeros(n) : nothing
    diag_parameter = parameterXML(id = diag_id, value = diag(X),
            lower = lower_bounds)

    n_offdiag = div(n * (n - 1), 2)
    offdiag_elements = zeros(n_offdiag)
    ind = 0
    for i = 1:n
        for j = (i + 1):n
            ind += 1
            offdiag_elements[ind] = X[i, j]
        end
    end

    offdiag_parameter = parameterXML(id = offdiag_id, value = offdiag_elements)

    return compoundSymmetricMatrixXML(diag_parameter, offdiag_parameter;
            kw_args...)
end

function compoundSymmetricMatrixXML(dim::Int, diag_id::String,
        offdiag_id::String; kw_args...)
    X = Diagonal(ones(dim))
    return compoundSymmetricMatrixXML(X, diag_id, offdiag_id; kw_args...)
end

function matrixParameterXML(parameters::Vector{<:AbstractGeneralizedXMLElement};
        id::StringOrNothing = nothing)
    return GeneralizedXMLElement("matrixParameter", id = id, children = parameters)
end

function matrixParameterXML(X::AbstractMatrix{Float64};
        id::StringOrNothing = nothing)

    n = size(X, 1)
    ids = isnothing(id) ? fill(nothing, n) : ["$id.$i" for i = 1:n]
    parameters = [parameterXML(id=ids[i], value=X[i, :]) for i = 1:n]
    return matrixParameterXML(parameters, id = id)
end

################################################################################
## mcmc
################################################################################

function mcmcXML(joint::GeneralizedXMLElement, operators::GeneralizedXMLElement,
        screen_log::GeneralizedXMLElement,
        file_log::GeneralizedXMLElement,
        chain_length::Int;
        auto_optimize::Bool = true)
    return GeneralizedXMLElement("mcmc", id = "mcmc",
        children = [joint, operators, screen_log, file_log],
        attributes = ["chainLength" => chain_length,
                "autoOptimize" => auto_optimize])
end

function mcmcXML(org::Organizer, operators::GeneralizedXMLElement;
            chain_length::Int,
            screen_logEvery::Int,
            file_logEvery::Int,
            file_name::String,
            overwrite::Bool = true)
    likelihood = GeneralizedXMLElement("likelihood", id = "likelihood",
            children = org.likelihoods)
    prior = GeneralizedXMLElement("prior", id = "prior", children = org.priors)
    joint = GeneralizedXMLElement("joint", id="joint",
            children = [likelihood, prior])

    default_loggables = [joint, likelihood, prior]
    screen_log = logXML(default_loggables, screen_logEvery)
    file_log = logXML([default_loggables; org.loggables], file_logEvery,
        file_name = file_name,
        overwrite = overwrite)

    return mcmcXML(joint, operators, screen_log, file_log, chain_length)
end

function logXML(elements::Vector{<:GeneralizedXMLElement}, log_every::Int;
        id::StringOrNothing = nothing,
        overwrite::Union{Nothing, Bool} = nothing,
        file_name::StringOrNothing = nothing)

    attrs = ["fileName" => file_name,
            "logEvery" => log_every,
            "overwrite" => overwrite]
    return GeneralizedXMLElement("log", id = id, children = elements,
            attributes = attrs)
end
