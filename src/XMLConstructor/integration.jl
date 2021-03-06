struct TextXMLElement <: MyXMLElement
    el::XMLElement
end

function name(el::TextXMLElement)
    return name(el.el)
end

function make_xml(el::TextXMLElement)
    # do nothing
end

function get_loggables(el::TextXMLElement)
    return [el.el]
end

function get_priors(el::TextXMLElement)
    return [el.el]
end


################################################################################
## Build from existing xml
################################################################################


function BEASTXMLElement(xml::String)
    xdoc = parse_file(xml)
    xroot = root(xdoc)
    if name(xroot) != bn.BEAST
        error("Not a BEAST xml.")
    end

    bx = BEASTXMLElement()

    for element in child_elements(xroot)
        el = parse_element(element)
        if !isnothing(el)
            add_child(bx, el)
        end
    end
    return bx
end

# const SPECIAL_PARSERS = Dict(bn.TAXA => parse_taxa)

function parse_element(el::XMLElement)
    nm = name(el)
    if nm == bn.TAXA
        return parse_taxa(el)
    elseif nm == "coalescentSimulator"
        return nothing
    elseif nm == bn.TREE_MODEL
        return nothing
    elseif nm == bn.OPERATORS
        return parse_operators(el)
    elseif nm == bn.MCMC
        return parse_mcmc(el)
    end
    return TextXMLElement(el)
end


struct EmptyDataXMLElement <: MyXMLElement
    taxa::Vector{String}
end

function parse_taxa(el::XMLElement)
    taxa_els = el[bn.TAXON]
    taxa = [get_id(x) for x in taxa_els]
    return EmptyDataXMLElement(taxa)
end

function parse_operators(el::XMLElement)
    ops = TextXMLElement[]
    for el in child_elements(el)
        push!(ops, TextXMLElement(el))
    end
    return OperatorsXMLElement(ops)
end

struct ParsedMCMCXMLElement <: MyXMLElement
    priors::Vector{TextXMLElement}
    likelihoods::Vector{TextXMLElement}
    file_loggables::Vector{TextXMLElement}
    tree_log::TextXMLElement
end


function parse_mcmc(el::XMLElement)
    j_el = find_element(el, "joint")
    prior_el = find_element(j_el, bn.PRIOR)
    priors = [TextXMLElement(cel) for cel in child_elements(prior_el)]

    like_el = find_element(j_el, bn.LIKELIHOOD)
    likes = [TextXMLElement(cel) for cel in collect(child_elements(like_el))]

    log_els = el[bn.LOG]
    file_loggables = nothing
    for lg in log_els
        id = get_id(lg)
        if id == "fileLog"
            c_els = child_elements(lg)
            file_loggables = TextXMLElement[]
            ind = 1
            for el in c_els
                if ind > 3
                    push!(file_loggables, TextXMLElement(el))
                end
                ind += 1
            end
        end
    end


    tree_log = TextXMLElement(find_element(el, "logTree"))
    j_el = find_element(tree_log.el, "joint")
    unlink(j_el)
    free(j_el)
    post_el = new_child(tree_log.el, bn.POSTERIOR)
    set_attribute(post_el, bn.IDREF, bn.POSTERIOR)

    return ParsedMCMCXMLElement(priors, likes, file_loggables, tree_log)
end








################################################################################
## Merge xml
################################################################################

# This will only work with some very specific circumnstances
#   1. Both have some el <: AbstractDataXMLElement

function findfirst_name(bx::BEASTXMLElement, nm::String)
    return findfirst(x -> name(x) == nm, bx.components)
end

function findall_name(bx::BEASTXMLElement, nm::String)
    return findall(x -> name(x) == nm, bx.components)
end

function merge_xml!(traitxml::BEASTXMLElement, seqxml::BEASTXMLElement)

    # merge taxa
    tx_trait = find_element(traitxml, DataXMLElement)
    # @show tx_trait
    tx_seq = find_element(seqxml, EmptyDataXMLElement)
    taxa_trait = tx_trait.taxa
    taxa_seq = tx_seq.taxa
    @assert length(taxa_trait) == length(taxa_seq)
    @assert Set(taxa_trait) == Set(taxa_seq) # TODO: more efficient?

    seq_start = findfirst(x -> typeof(x) <: EmptyDataXMLElement, seqxml.components) + 1

    seq_ops_ind = findfirst(x -> typeof(x) <: OperatorsXMLElement, seqxml.components)
    seq_stop = seq_ops_ind - 1

    trait_ind = findfirst(x -> typeof(x) <: TreeModelXMLElement, traitxml.components) + 1

    for i = seq_start:seq_stop
        add_child(traitxml, seqxml.components[i], trait_ind)
        trait_ind += 1
    end

    trait_ops_ind = findfirst(x -> typeof(x) <: OperatorsXMLElement, traitxml.components)
    merge_operators!(traitxml.components[trait_ops_ind],
                     seqxml.components[seq_ops_ind])

    seqmc = find_element(seqxml, ParsedMCMCXMLElement)
    traitmc = find_element(traitxml, MCMCXMLElement)
    merge_mcmc!(traitmc, seqmc)
    add_loggables(traitxml, LoggablesXMLElement(seqmc.file_loggables))

    newick_el = find_element(traitxml, NewickXMLElement)
    newick_el.fix_tree = false


end