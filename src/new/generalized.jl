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
        add_text(el, content)
    end

    xml.defined = true

    return el
end

function find_elements(xml::GeneralizedXMLElement;
        name::String = "",
        attributes::Dict{String, String} = Dict{String, String}())


    check_name = !isempty(name)
    check_attributes = !isempty(attributes)

    if !(check_name || check_attributes)
        error("must supply either name or attributes (or both) to match")
    end

    matches = GeneralizedXMLElement[]



    for child in xml.children
        if check_name
            if child.name != name
                continue
            end
        end

        if check_attributes
            if !matches_attributes(child, attributes)
                continue
            end
        end
        push!(matches, child)
    end

    return matches
end

function matches_attributes(xml::GeneralizedXMLElement,
        attributes::Dict{String, String})
    valid = true
    for (k, v) in attributes
        has_attribute = false
        for (xk, xv) in xml.attributes
            if xk == k
                has_attribute = true
                if xv != v
                    valid = false
                    break
                end
            end
            if !valid
                break
            end
        end
        if !has_attribute
            valid = false
            break
        end
    end

    return valid
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



