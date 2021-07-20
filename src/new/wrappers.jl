struct Organizer
    elements::Vector{<:GeneralizedXMLElement}
    likelihoods::Vector{<:GeneralizedXMLElement}
    priors::Vector{<:GeneralizedXMLElement}
    loggables::Vector{<:GeneralizedXMLElement}
end

function Organizer(elements::Vector{<:GeneralizedXMLElement};
        likelihoods::Vector{<:GeneralizedXMLElement} = GeneralizedXMLElement[],
        priors::Vector{<:GeneralizedXMLElement} = GeneralizedXMLElement[],
        loggables::Vector{<:GeneralizedXMLElement} = GeneralizedXMLElement[])

    return Organizer(elements, likelihoods, priors, loggables)
end

const CONCAT_FIELDS = [:elements, :likelihoods, :priors, :loggables]

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

