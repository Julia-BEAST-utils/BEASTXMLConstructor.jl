mutable struct GeneralDataType <: MyXMLElement
    el::XMLOrNothing
    states::Vector{Int}
    id::String

    function GeneralDataType(states::Vector{Int}, id::String)
        return new(nothing, states, id::String)
    end
end

function GeneralDataType(n::Int; id::String = "discreteStates")
    return GeneralDataType(collect(0:(n - 1)), id)
end

function name(::GeneralDataType)
    return bn.GENERAL_DATA_TYPE
end

function make_xml(gdt::GeneralDataType)
    el = new_element(name(gdt))
    set_attribute(el, bn.ID, gdt.id)

    for i = 1:length(gdt.states)
        state_el = new_child(el, bn.STATE)
        set_attribute(state_el, bn.CODE, gdt.states[i])
    end

    gdt.el = el
    return el
end


################################################################################
## aligments
################################################################################

const MISSING_SEQUENCE = 9

mutable struct Alignment <: MyXMLElement
    el::XMLOrNothing
    data_type::MyXMLElement
    taxa::Vector{<:AbstractString}
    sequences::Matrix{Int}
    id::String

    function Alignment(data_type::MyXMLElement,
                        taxa::Vector{<:AbstractString},
                        sequences::Matrix{Int};
                        id::String = bn.ALIGNMENT)
        @assert length(taxa) == size(sequences, 1)
        return new(nothing, data_type, taxa, sequences, id)
    end
end

function name(::Alignment)
    return bn.ALIGNMENT
end

function make_xml(al::Alignment)
    el = new_element(name(al))
    add_ref_el(el, al.data_type)
    for i = 1:size(al.sequences, 1)
        seq_el = new_child(el, bn.SEQUENCE)
        taxon_el = new_child(seq_el, bn.TAXON)
        set_attribute(taxon_el, bn.IDREF, al.taxa[i])
        seq_text = join(al.sequences[i, :], ' ')
        add_text(seq_el, seq_text)
    end

    al.el = el
    return el
end

function get_num_classes(alignment::Alignment)
    @unpack sequences = alignment
    p = size(alignment.sequences, 2)
    return length.([setdiff(unique(sequences[:, i]), MISSING_SEQUENCE) for i = 1:p])
end

################################################################################
## patterns
################################################################################

mutable struct Patterns <: MyXMLElement
    el::XMLOrNothing
    alignment::Alignment
    id::String

    function Patterns(alignment::Alignment; id::String = bn.PATTERNS)
        return new(nothing, alignment, id)
    end
end

function name(::Patterns)
    return bn.PATTERNS
end

function make_xml(patterns::Patterns)
    @unpack alignment, id = patterns

    el = new_element(name(patterns))
    set_attribute(el, bn.ID, id)
    set_attribute(el, bn.FROM, "1")
    set_attribute(el, bn.UNIQUE, bn.FALSE)
    add_ref_el(el, alignment)
    patterns.el = el

    return el
end

function get_num_classes(patterns::Patterns)
    return get_num_classes(patterns.alignment)
end

################################################################################
## likelihood
################################################################################

mutable struct OrderedLatentLiability <: MyXMLElement
    el::XMLOrNothing
    patterns::Patterns
    treeModel::TreeModelXMLElement
    threshold::CompoundParameter
    id::String
    trait_name::String

    function OrderedLatentLiability(patterns::Patterns,
                                    treeModel::TreeModelXMLElement,
                                    threshold::CompoundParameter;
                                    id::String = bn.ORDERED_LATENT_LIABILITY,
                                    trait_name::String = treeModel.param_names[1])

        return new(nothing, patterns, treeModel, threshold, id, trait_name)
    end
end

function OrderedLatentLiability(patterns::Patterns,
                                treeModel::TreeModelXMLElement;
                                id::String = bn.ORDERED_LATENT_LIABILITY,
                                trait_name::String = treeModel.param_names[1])
    n_states = get_num_classes(patterns)
    need_thresh = findall(x -> x > 2, n_states)

    params = Parameter[]
    if length(need_thresh) == 0
        push!(params, Parameter([0.0]))
    else
        for ind in need_thresh
            p_id = "$trait_name$ind.threshold"
            thresholds = default_thresholds(n_states[ind])
            push!(params, Parameter(thresholds, p_id))
        end
    end

    cp = CompoundParameter(params, "$trait_name.threshold")
    return OrderedLatentLiability(patterns, treeModel, cp, id = id, trait_name = trait_name)

end

function name(::OrderedLatentLiability)
    return bn.ORDERED_LATENT_LIABILITY
end

function make_xml(oll::OrderedLatentLiability)
    @unpack patterns, treeModel, threshold, id = oll
    el = new_element(name(oll))
    add_ref_el(el, patterns)
    add_ref_el(el, treeModel)
    tip_el = new_child(el, bn.TIP_TRAIT)
    add_ref_el(tip_el, Parameter(Float64[], oll.trait_name))

    thresh_el = new_child(el, bn.THRESHOLD)
    add_child(thresh_el, make_xml(threshold))

    classes_el = new_element(bn.NUM_CLASSES)
    classes = get_num_classes(patterns)
    class_param = Parameter(Float64.(classes))
    add_child(classes_el, make_xml(class_param))

    oll.el = el
    return el
end


################################################################################
## operator
################################################################################

mutable struct ExtendedLatendLiabilityOperator <: OperatorXMLElement
    el::XMLOrNothing
    trait_likelihood::TraitLikelihoodXMLElement
    latent_liability::OrderedLatentLiability
    weight::Float64

    function ExtendedLatendLiabilityOperator(
            trait_likelihood::TraitLikelihoodXMLElement,
            latent_liability::OrderedLatentLiability;
            weight::Float64 = 1.0)

        return new(nothing, trait_likelihood, latent_liability, weight)
    end

end

function name(::ExtendedLatendLiabilityOperator)
    return bn.EXTENDED_GIBBS_OPERATOR
end

function make_xml(elo::ExtendedLatendLiabilityOperator)
    el = new_element(name(elo))

    set_attribute(el, bn.WEIGHT, elo.weight)
    add_ref_el(el, elo.trait_likelihood)
    add_ref_el(el, elo.latent_liability)
    elo.el = el
    return el
end

################################################################################
## helpers
################################################################################

function default_thresholds(states::Int)
    return convert.(Float64, collect(0:(states - 2)))
end

function find_state(x::Float64, thresholds::Vector{Float64})
    n = length(thresholds)

    if isnan(x)
        if n >= MISSING_SEQUENCE
            error("Cannot accomodate more than $MISSING_SEQUENCE patterns current missing code.")
        end
        return MISSING_SEQUENCE
    end

    state = 0
    for i = 1:n
        if x < thresholds[i]
            return i - 1
        end
    end
    return n
end

function data_to_alignment(data::Matrix{Float64}, states::Vector{Int};
            thresholds::Vector{Vector{Float64}} = default_thresholds.(states),
            bin_discrete::Bool = false)

    n, p = size(data)
    sequences = zeros(Int, n, p)

    for j = 1:p
        n_states = states[j]
        if n_states > 1
            if bin_discrete
                for i = 1:n
                    sequences[i, j] = find_state(data[i, j], thresholds[j])
                end
            else
                sequences[:, j] .= Int.(data[:, j])
            end
        end
    end

    return sequences
end

################################################################################
## latent liability
################################################################################

function add_latent_liability(bx::BEASTXMLElement, trait_names::Vector{String}, states::Vector{Int})
    m = length(trait_names)
    gdt = GeneralDataType(maximum(states))
    ops_ind = find_element_ind(bx, OperatorsXMLElement)

    add_child(bx, gdt, ops_ind)

    for nm in trait_names
        latent_liablity = latent_liablity_elements(bx, nm, states, gdt)
    end


end

function latent_liablity_elements(bx::BEASTXMLElement, trait_name::String, states::Vector{Int}, gdt::GeneralDataType)

    data_el = find_element(bx, DataXMLElement)
    treeModel = find_element(bx, TreeModelXMLElement)
    trait_likelihood = find_element(bx, TraitLikelihoodXMLElement)

    taxa = data_el.taxa
    trait_ind = findfirst(isequal(trait_name), data_el.trait_names)
    data = data_el.data_mats[trait_ind]

    sequence = data_to_alignment(data, states)
    alignment = Alignment(gdt, taxa, sequence, id=trait_name * ".alignment")
    pattern = Patterns(alignment, id=trait_name * ".patterns")

    param_name = treeModel.param_names[findfirst(isequal(trait_name), treeModel.node_traits)]

    likelihood = OrderedLatentLiability(pattern, treeModel, id = trait_name * ".latentLiability", trait_name = param_name)

    operator = ExtendedLatendLiabilityOperator(trait_likelihood, likelihood)

    ops_ind = find_element_ind(bx, OperatorsXMLElement)
    operators = bx.components[ops_ind]
    for el in [alignment, pattern, likelihood]
        add_child(bx, el, ops_ind)
        ops_ind += 1
    end

    push!(operators.els, operator)



    return likelihood
end

