mutable struct GeneralDataType <: MyXMLElement
    el::XMLOrNothing
    states::Vector{Char}
    id::String

    function GeneralDataType(states::Vector{Char}, id::String)
        return new(nothing, states, id::String)
    end
end

# function GeneralDataType(n::Int; id::String = "discreteStates")
#     return GeneralDataType(collect(0:(n - 1)), id)
# end

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
    sequences::Matrix{Char}
    id::String

    function Alignment(data_type::MyXMLElement,
                        taxa::Vector{<:AbstractString},
                        sequences::Matrix{Char};
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
    set_attribute(el, bn.ID, al.id)
    add_ref_el(el, al.data_type)
    for i = 1:size(al.sequences, 1)
        seq_el = new_child(el, bn.SEQUENCE)
        taxon_el = new_child(seq_el, bn.TAXON)
        set_attribute(taxon_el, bn.IDREF, al.taxa[i])
        seq_brackets = ["$x" for x in al.sequences[i, :]]
        seq_text = join(seq_brackets, ' ')
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
            push!(params, Parameter(thresholds, p_id, lower=0.0))
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

    classes_el = new_child(el, bn.NUM_CLASSES)
    classes = get_num_classes(patterns)
    class_param = Parameter(Float64.(classes))
    add_child(classes_el, make_xml(class_param))

    set_attribute(el, bn.N_TRAITS, string(length(classes)))
    set_attribute(el, bn.N_DATA, "1")
    set_attribute(el, bn.ID, oll.id)


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
    return convert.(Float64, collect(1:(states - 2)))
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

function data_to_alignment(data::Matrix{Float64}, states::Vector{Int},
            states_dict::Dict{Int, Char};
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
                for i = 1:n
                    if isnan(data[i, j])
                        sequences[i, j] = MISSING_SEQUENCE
                    else
                        sequences[i, j] = Int(data[i, j])
                    end
                end
            end
        end
    end

    str_sequences = [states_dict[sequences[i, j]] for i = 1:n, j = 1:p]

    return str_sequences

end

################################################################################
## latent liability
################################################################################


const ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"

function add_latent_liability(bx::BEASTXMLElement, trait_names::Vector{String}, states::Vector{Int})
    m = length(trait_names)
    mx = maximum(states)
    if mx >= length(ALPHABET)
        error("Cannot accommodated more than $(length(ALPHABET)) unique discrete states.")
    end

    str_states = [ALPHABET[i] for i = 1:(mx + 1)]  # need to add 1 for zero state
    state_dict = Dict(i => str_states[i + 1] for i = 0:mx)

    gdt = GeneralDataType(str_states, "discreteStates")
    ops_ind = find_element_ind(bx, OperatorsXMLElement)

    add_child(bx, gdt, ops_ind)

    for nm in trait_names
        latent_liablity = latent_liablity_elements(bx, nm, states, state_dict,
                gdt)
    end
end

function add_latent_liability(bx::BEASTXMLElement, discrete_inds::Vector{Int})

    if length(discrete_inds) == 0
        return nothing
    end

    data_el = find_element(bx, DataXMLElement)
    @unpack data_mats = data_el
    ps = [size(x, 2) for x in data_mats]
    @assert length(unique(ps)) == 1
    p = ps[1]
    states = ones(Int, p)

    for data in data_el.data_mats
        for ind in discrete_inds
            ind_states = setdiff(unique(@view data[:, ind]), NaN)
            states[ind] = max(Int(maximum(ind_states)) + 1, states[ind])
        end
    end

    set_factor_precision_indices(bx, setdiff(1:p, discrete_inds))



    add_latent_liability(bx, data_el.trait_names, states)
end



function latent_liablity_elements(bx::BEASTXMLElement, trait_name::String,
                                  states::Vector{Int},
                                  states_dict::Dict{Int, Char},
                                  gdt::GeneralDataType)

    data_el = find_element(bx, DataXMLElement)
    treeModel = find_element(bx, TreeModelXMLElement)
    trait_likelihood = find_element(bx, TraitLikelihoodXMLElement)
    int_fac = trait_likelihood.extension_el

    if BEASTXMLConstructor.get_trait_name(int_fac) != trait_name
        flpd = find_element(bx, FactorLogPredictiveDensity)
        trait_likelihood = flpd.like
    end

    @assert BEASTXMLConstructor.get_trait_name(trait_likelihood.extension_el) == trait_name

    precs = get_precision(int_fac)
    set_value(precs, fill(2.0, length(precs)))

    taxa = data_el.taxa
    trait_ind = findfirst(isequal(trait_name), data_el.trait_names)
    data = data_el.data_mats[trait_ind]

    sequence = data_to_alignment(data, states, states_dict)
    alignment = Alignment(gdt, taxa, sequence, id=trait_name * ".alignment")
    pattern = Patterns(alignment, id=trait_name * ".patterns")

    param_name = treeModel.param_names[findfirst(isequal(trait_name), treeModel.node_traits)]

    likelihood = OrderedLatentLiability(pattern, treeModel, id = trait_name * ".latentLiability", trait_name = param_name)

    operator = ExtendedLatendLiabilityOperator(trait_likelihood, likelihood)
    thresh_ops = MyXMLElement[]
    if (maximum(get_num_classes(likelihood.patterns)) > 2)
        for param in likelihood.threshold.parameters
            push!(thresh_ops, ScaleOperator(param))
        end
    end
    # threshold_operator = ScaleOperator(likelihood.threshold)

    ops_ind = find_element_ind(bx, OperatorsXMLElement)
    operators = bx.components[ops_ind]
    for el in [alignment, pattern, likelihood]
        add_child(bx, el, ops_ind)
        ops_ind += 1
    end

    push!(operators.els, operator)
    for op in thresh_ops
        push!(operators.els, op)
    end

    mcmc_el = find_element(bx, MCMCXMLElement)
    BEASTXMLConstructor.add_likelihood!(mcmc_el, likelihood)
    BEASTXMLConstructor.add_loggable(bx, likelihood.threshold, already_made = true)
    BEASTXMLConstructor.add_loggable(bx, Parameter(Float64[], param_name), already_made = true)



    return likelihood
end

