mutable struct TreeModelXMLElement <: MyXMLElement
    el::XMLOrNothing
    node_traits::Vector{String}
    trait_dims::Vector{Int}
    param_names::Vector{String}
    newick_el::NewickXMLElement
    fixed_tree::Bool
    id::String

    function TreeModelXMLElement(newick_el::NewickXMLElement,
                dim::Int;
                trait::String = bn.DEFAULT_TRAIT_NAME)
        param_name = bn.LEAF_TRAITS
        return new(nothing, [trait], [dim], [param_name], newick_el, true,
                    bn.TREE_MODEL)
    end

    function TreeModelXMLElement(newick_el::NewickXMLElement,
                node_traits::Vector{String},
                trait_dims::Vector{Int},
                param_names::Vector{String}
                )
        return new(nothing, node_traits, trait_dims, param_names, newick_el, true)
    end

    function TreeModelXMLElement(newick_el::NewickXMLElement,
                                 node_traits::Vector{String},
                                 trait_dims::Vector{Int}
                                )
        param_names = ["$trait.leafTrait" for trait in node_traits]
        return TreeModelXMLElement(newick_el, node_traits, trait_dims, param_names)
    end
end

function name(::TreeModelXMLElement)
    return bn.TREE_MODEL
end

function get_id(tm::TreeModelXMLElement)
    return tm.id
end

function make_xml(tl::TreeModelXMLElement)
    make_xml(tl.newick_el)
    tl.el =  make_treeModel(tl.newick_el.el,
                    tl.node_traits,
                    tl.trait_dims,
                    tl.param_names,
                    id=tl.id)
end

function make_treeModel(newick_el::XMLElement,
        trait_names::Vector{String},
        trait_dims::Vector{Int},
        param_names::Vector{String};
        fixed_tree::Bool = false,
        id::String = bn.TREE_MODEL)


    @assert name(newick_el) == bn.NEWICK

    el = new_element(bn.TREE_MODEL)
    set_attribute(el, bn.ID, id)
    set_attribute(el, bn.FIX_HEIGHTS, bn.TRUE)

    if fixed_tree
        set_attribute(el, bn.FIXED_TREE, bn.TRUE)
    end

    tree_el = new_element(bn.TREE)
    newick_id = attribute(newick_el, bn.ID)
    set_attribute(tree_el, bn.IDREF, newick_id)
    add_child(el, tree_el)

    rh_el = new_element(bn.ROOT_HEIGHT)
    add_parameter_id(rh_el, "$(bn.TREE_MODEL).$(bn.ROOT_HEIGHT)")
    add_child(el, rh_el)

    add_nodeHeights(el, "$(bn.TREE_MODEL).$(bn.INTERNAL_NODE_HEIGHTS)",
            internal = true)
    add_nodeHeights(el, "$(bn.TREE_MODEL).$(bn.ALL_INTENAL_NODE_HEIGHTS)",
            internal = true, root = true)

    for i = 1:length(trait_names)
        nt_el = new_element(bn.NODE_TRAITS)
        attrs = [(bn.ROOT_NODE, bn.FALSE),
                (bn.INTERNAL_NODES, bn.FALSE),
                (bn.LEAF_NODES, bn.TRUE),
                (bn.AS_MATRIX, bn.TRUE),
                (bn.TRAIT_DIM, string(trait_dims[i])),
                (bn.NAME, trait_names[i])]
        set_attributes(nt_el, attrs)
        add_parameter_id(nt_el, param_names[i])
        add_child(el, nt_el)
    end
    return el
end

function add_nodeHeights(pel::XMLElement, param_id::String;
        internal::Bool = false, root::Bool = false)

    el = new_element(bn.NODE_HEIGHTS)
    if internal
        set_attribute(el, bn.INTERNAL_NODES, bn.TRUE)
    end
    if root
        set_attribute(el, bn.ROOT_NODE, bn.TRUE)
    end
    add_parameter_id(el, param_id)
    add_child(pel, el)
    return el
end

function add_leaf_param!(txl::TreeModelXMLElement, trait_name::String, trait_dim::Int)
    push!(txl.node_traits, trait_name)
    push!(txl.trait_dims, trait_dim)
    push!(txl.param_names, "$(trait_name).leaf")
end
