struct Organizer
    elements::Vector{<:GeneralizedXMLElement}
    likelihoods::Vector{<:GeneralizedXMLElement}
    priors::Vector{<:GeneralizedXMLElement}
    loggables::Vector{<:GeneralizedXMLElement}
    partial_providers::Vector{<:GeneralizedXMLElement}
end

function Organizer(elements::Vector{<:GeneralizedXMLElement};
        likelihoods::Vector{<:GeneralizedXMLElement} = GeneralizedXMLElement[],
        priors::Vector{<:GeneralizedXMLElement} = GeneralizedXMLElement[],
        loggables::Vector{<:GeneralizedXMLElement} = GeneralizedXMLElement[],
        partial_providers::Vector{<:GeneralizedXMLElement} = GeneralizedXMLElement[])

    return Organizer(elements, likelihoods, priors, loggables, partial_providers)
end

const CONCAT_FIELDS = fieldnames(Organizer)

import Base.vcat
function vcat(orgs::Organizer...)
    n_fields = length(CONCAT_FIELDS)
    new_fields = [GeneralizedXMLElement[] for i = 1:n_fields]
    for i = 1:length(orgs)
        for j = 1:n_fields
            new_fields[j] = [new_fields[j]; getfield(orgs[i], CONCAT_FIELDS[j])]
        end
    end
    return Organizer(new_fields...) #TODO: will need to modify for more complex struct
end

import Base.push!
function push!(org::Organizer, el::GeneralizedXMLElement;
        loggable::Bool = false, likelihood::Bool = false, prior::Bool = false)
    push!(org.elements, el)
    if loggable
        push!(org.loggables, el)
    end
    if likelihood
        push!(org.likelihoods, el)
    end
    if prior
        push!(org.priors, el)
    end
    return org
end

function find_element(org::Organizer, name::String)
    i = findfirst(x -> x.name == name, org.elements)
    if isnothing(i)
        return nothing
    end
    return org.elements[i]
end

