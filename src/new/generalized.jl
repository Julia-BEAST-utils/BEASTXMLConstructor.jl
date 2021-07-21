abstract type AbstractGeneralizedXMLElement end
const StringOrNothing = Union{String, Nothing}

################################################################################
## Utils
################################################################################
function format_string(x)
    return string(x)
end

function format_string(x::AbstractVector{<:Real})
    return join(x, ' ')
end

function format_string(X::AbstractMatrix{<:Real})
    return join([join(X[i, :], ' ') for i = 1:size(X, 1)], '\n')
end

function format_string(::Nothing)
    return nothing
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



function matches_attributes(::PassthroughXMLElement, ::Dict{String, String})
    return false
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
    formatted_attributes = Pair{String, String}[]
    for i = 1:p
        if !isnothing(attributes[i][2])
            push!(formatted_attributes, attributes[i][1] =>
                    format_string(attributes[i][2]))
        end
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
        @show xml.name
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
        add_text(el, content)
    end

    xml.defined = true

    return el
end

function find_elements(xml::GeneralizedXMLElement;
        names::Vector{String} = String[],
        name::String = "",
        attributes::Dict{String, String} = Dict{String, String}(),
        passthrough::Bool = false)


    check_name = !isempty(name) || !isempty(names)
    if check_name
        if !isempty(name) && !(isempty(names))
            error("only supply on the 'name' or 'names' keyword arguments")
        elseif isempty(names)
            names = [name]
        end
    end
    check_attributes = !isempty(attributes)

    if !(check_name || check_attributes)
        error("must supply either name or attributes (or both) to match")
    end

    matches = GeneralizedXMLElement[]



    for child in xml.children
        if element_matches(child, names, attributes,
                check_name = check_name, check_attributes = check_attributes,
                passthrough = passthrough)
            add_match!(matches, child)
        end
    end

    return matches
end

function add_match!(matches::Vector{GeneralizedXMLElement}, match::GeneralizedXMLElement)
    push!(matches, match)
end

function add_match!(matches::Vector{GeneralizedXMLElement}, match::PassthroughXMLElement)
    if length(match.children) != 1
        error("can only currenlty add passthrough elements with one child")
    end
    push!(matches, match.children[1])
end

# function find_parent_elements(xml::GeneralizedXMLElement, name::String)
#     matches = GeneralizedXMLElement[]
#     for child in children(xml)



function element_matches(xml::GeneralizedXMLElement, names::Vector{String},
            attrs::Dict{String, String};
            check_name::Bool,
            check_attributes::Bool,
            passthrough::Bool = false)

    if passthrough
        return false
    end

    @assert check_name || check_attributes

    if check_name && !(xml.name in names)
        return false
    end
    if check_attributes && !matches_attributes(xml, attrs)
        return false
    end

    return true
end

function element_matches(xml::PassthroughXMLElement, names::Vector{String},
        attrs::Dict{String, String};
        check_name::Bool,
        check_attributes::Bool,
        passthrough::Bool = false)

    if !passthrough
        return false
    end
    if check_attributes
        error("cannot currently check attributes for PassthroughXMLElement")
    end
    @assert check_name

    if !(xml.name in names)
        return false
    end

    return true
end


function matches_attributes(xml::GeneralizedXMLElement,
        attributes::Dict{String, String})
    for (k, v) in attributes
        has_attribute = false
        for (xk, xv) in xml.attributes
            if xk == k
                has_attribute = true
                if xv != v
                    return false
                end
            end
        end
        if !has_attribute
            return false
        end
    end

    return true
end


function find_element(xml::GeneralizedXMLElement; kw_args...)
    matches = find_elements(xml; kw_args...)
    if length(matches) == 0
        error("no elements found matching criteria")
    elseif length(matches) > 1
        error("more than one element found matching creteria")
    end

    return matches[1]
end

function get_attribute(xml::GeneralizedXMLElement, attribute::String)
    ind = findfirst(x -> x[1] == attribute, xml.attributes)
    if isnothing(ind)
        error("attribute '$attribute' not found in $(xml.name) element")
    end
    return xml.attributes[ind][2]
end

function get_id(xml::GeneralizedXMLElement)::String
    if isnothing(xml.id)
        error("element $(xml.name) has no id")
    end
    return xml.id
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



