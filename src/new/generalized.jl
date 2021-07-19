abstract type AbstractGeneralizedXMLElement end
const StringOrNothing = Union{String, Nothing}

################################################################################
## Utils
################################################################################
function format_string(x)
    return string(x)
end

function format_string(x::AbstractVector{Float64})
    return join(x, ' ')
end

function format_string(X::AbstractMatrix{Float64})
    return join([join(X[i, :], ' ') for i = 1:size(X, 1)], '\n')
end

function format_string(::Nothing)
    return nothing
end

################################################################################
## GeneralizedXMLElement
################################################################################

mutable struct GeneralizedXMLElement{T<: StringOrNothing, S<:StringOrNothing} <: AbstractGeneralizedXMLElement
    name::String
    id::T
    attributes::Array{Pair{String, String}}
    children::Vector{<:AbstractGeneralizedXMLElement}
    content::S
    defined::Bool
end

function GeneralizedXMLElement(name::String, id::StringOrNothing,
                               attributes::Array{Pair{String, String}},
                               children::Vector{<:AbstractGeneralizedXMLElement},
                               content::StringOrNothing)

    GeneralizedXMLElement(name, id, attributes, children, content, false)
end

function GeneralizedXMLElement(name::String, id::StringOrNothing,
        attributes::Array{<:Pair{String, <:Any}},
        children::Vector{<:AbstractGeneralizedXMLElement},
        content::Any)

    p = length(attributes)
    formatted_attributes = Vector{Pair{String, String}}(undef, p)
    for i = 1:p
        formatted_attributes[i] = attributes[i][1] => format_string(attributes[i][2])
    end

    formatted_content = format_string(content)

    GeneralizedXMLElement(name, id, formatted_attributes, children,
        formatted_content, false)
end


function GeneralizedXMLElement(name::String;
        id::StringOrNothing = nothing,
        attributes::Array{<:Pair{String, <:Any}} = Pair{String, String}[],
        children::Vector{<:AbstractGeneralizedXMLElement} = GeneralizedXMLElement[],
        child::Union{Nothing, AbstractGeneralizedXMLElement} = nothing,
        content::Any = nothing)

    if !isnothing(child)
        if length(children) > 0
            error("cannot supply both 'child' and 'children' keyword arguments")
        else
            children = [child]
        end
    end

    return GeneralizedXMLElement(name, id, attributes, children, content)
end



function already_defined(xml::GeneralizedXMLElement)
    return xml.defined
end

function has_id(::GeneralizedXMLElement{Nothing, T}) where T
    return false
end

function has_id(::GeneralizedXMLElement{String, T}) where T
    return true
end


function reference_element(xml::GeneralizedXMLElement)
    if !already_defined(xml)
        error("Cannot reference element that has not alreday been constructed.")
    elseif !has_id(xml)
        error("Cannot reference element without id attribute.")
    end

    el = new_element(xml.name)
    set_attribute(el, "idref", xml.id)
    return el
end

function make_element(xml::GeneralizedXMLElement)::XMLElement
    @unpack name, id, attributes, children, content, defined = xml
    el = new_element(name)
    if defined
        return reference_element(xml)
    end

    if has_id(xml)
        set_attribute(el, "id", id)
    end

    set_attributes(el, attributes)
    for child in children
        add_child(el, make_element(child))
    end

    if !isnothing(content)
        @show content
        add_text(el, content)
    end

    xml.defined = true

    return el
end

################################################################################
## PassthroughXMLElement
################################################################################

struct PassthroughXMLElement <: AbstractGeneralizedXMLElement
    name::String
    children::Vector{<:AbstractGeneralizedXMLElement}
end

function PassthroughXMLElement(name::String, child::AbstractGeneralizedXMLElement)
    return PassthroughXMLElement(name, [child])
end

function make_element(xml::PassthroughXMLElement)
    el = new_element(xml.name)
    for child in xml.children
        add_child(el, make_element(child))
    end
    return el
end

################################################################################
## XMLDocument
################################################################################

struct BEASTXMLDocument
    elements::Vector{<:AbstractGeneralizedXMLElement}
end

function make_xml(document::BEASTXMLDocument)
    #TODO: reset all 'defined' fields
    bxml = GeneralizedXMLElement("beast", children = document.elements)
    beast_element = make_element(bxml)

    xdoc = XMLDocument()
    set_root(xdoc, beast_element)
    return xdoc
end

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

const DataPairs = Vector{Pair{String, Matrix{Float64}}}

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

function treeModelXML(tree::GeneralizedXMLElement, data::DataPairs)
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

    return GeneralizedXMLElement("treeModel", children = children)
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
    return GeneralizedContinuousTraitModel(id = id, children = parameters)
end

function matrixParameterXML(X::AbstractMatrix{Float64};
        id::StringOrNothing = nothing)

    n = size(X, 1)
    ids = isnothing(id) ? fill(nothing, n) : ["$id.$i" for i = 1:n]
    parameters = [parameterXML(id=ids[i], value=X[i, :]) for i = 1:n]
    return matrixParameterXML(parameters, id = id)
end




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
