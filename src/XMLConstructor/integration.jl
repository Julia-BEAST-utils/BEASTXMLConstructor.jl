struct TextXMLElement <: MyXMLElement
    el::XMLElement
end

function name(el::TextXMLElement)
    return name(el.el)
end

function make_xml(::TextXMLElement)
    # do nothing
end

function get_loggables(el::TextXMLElement)
    return [el.el]
end

function get_priors(el::TextXMLElement)
    return [el]
end

function add_ref_el(el::XMLElement, txt::TextXMLElement)
    add_child(el, txt.el)
end

function add_ref_els(el::XMLElement, txts::Vector{TextXMLElement})
    for txt in txts
        add_ref_el(el, txt)
    end
end

function get_id(el::TextXMLElement)
    id = attribute(el.el, bn.ID)
    if isnothing(id)
        id = attribute(el.el, bn.IDREF)
    end
    if isnothing(id)
        error("no id or idref")
    end

    return id
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
    taxa_parsed = false
    for element in child_elements(xroot)
        el = parse_element(element, taxa_parsed = taxa_parsed)
        if typeof(el) <: EmptyDataXMLElement
            taxa_parsed = true
        end
        if !isnothing(el)
            add_child(bx, el)
        end
    end
    return bx
end

# const SPECIAL_PARSERS = Dict(bn.TAXA => parse_taxa)

const DO_NOT_PARSE = ["coalescentSimulator", bn.TREE_MODEL, "rescaledTree", bn.NEWICK, "newick"]

function parse_element(el::XMLElement; taxa_parsed = false)
    nm = name(el)
    if nm == bn.TAXA && !taxa_parsed
        return parse_taxa(el)
    elseif nm in DO_NOT_PARSE
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
    file_log::TextXMLElement
    tree_log::TextXMLElement
    file_loggables::Vector{TextXMLElement}
end


function parse_mcmc(el::XMLElement)
    j_el = find_element(el, "joint")
    prior_el = find_element(j_el, bn.PRIOR)
    priors = [TextXMLElement(cel) for cel in child_elements(prior_el)]

    like_el = find_element(j_el, bn.LIKELIHOOD)
    likes = [TextXMLElement(cel) for cel in collect(child_elements(like_el))]

    file_log = nothing
    log_els = el[bn.LOG]
    file_loggables = nothing
    for lg in log_els
        id = get_id(lg)
        if id == "fileLog"
            set_attribute(lg, bn.ID, "fileLogSequence")
            file_log = TextXMLElement(lg)
            c_els = child_elements(lg)
            file_loggables = TextXMLElement[]
            ind = 1
            for el in c_els
                if ind > 3 #first thre are alwasy posterior, likelihood, prior
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

    return ParsedMCMCXMLElement(priors, likes, file_log, tree_log, file_loggables)
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

function merge_xml!(traitxml::BEASTXMLElement, seqxml::BEASTXMLElement; separate_logs::Bool = true)

    # merge taxa
    tx_trait = find_element(traitxml, DataXMLElement)
    # @show tx_trait
    tx_seq = find_element(seqxml, EmptyDataXMLElement)
    taxa_trait = tx_trait.taxa
    taxa_seq = tx_seq.taxa
    if Set(taxa_trait) != Set(taxa_seq)
        only_trait = setdiff(taxa_trait, taxa_seq)
        only_seq = setdiff(taxa_seq, taxa_trait)
        if length(only_seq) > 0
            error("The following taxa are present in the sequence xml but not the data:\n\t" *
                  join(only_seq, ' '))
        else
            error("The following taxa are present in the data but not the sequence xml:\n\t" *
                  join(only_trait, ' '))
        end
    end

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
    merge_mcmc!(traitmc, seqmc, separate_logs = separate_logs)

    # display(traitxml)
    # add_loggables(traitxml, LoggablesXMLElement(seqmc.file_loggables))
    # display(traitxml)
    # error()

    newick_el = find_element(traitxml, NewickXMLElement)
    newick_el.fix_tree = false
end