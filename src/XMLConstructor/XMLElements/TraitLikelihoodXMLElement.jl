mutable struct TraitLikelihoodXMLElement <: MyXMLElement
    el::XMLOrNothing
    mbd_el::MBDXMLElement
    treeModel_el::TreeModelXMLElement
    extension_el::Union{ModelExtensionXMLElement, Nothing}
    μ::Vector{Float64}
    pss::Float64
    attrs::Dict{String, String}
    xml_name::String


end

function TraitLikelihoodXMLElement(mbd_el::MBDXMLElement,
                    treeModel_el::TreeModelXMLElement,
                    extension_el::Union{ModelExtensionXMLElement, Nothing})
    μ = zeros(size(mbd_el.precision, 1))
    pss = 0.001
    trait_name = treeModel_el.node_traits[1]
    attrs = deepcopy(DEFAULT_TRAITLIKELIHOOD_ATTRS)
    attrs[bn.TRAIT_NAME] = trait_name
    return TraitLikelihoodXMLElement(nothing, mbd_el, treeModel_el,
                extension_el, μ, pss, attrs,
                bn.TRAIT_DATA_LIKELIHOOD)
end

function name(tl::TraitLikelihoodXMLElement)
    return tl.xml_name
end

function get_id(tl::TraitLikelihoodXMLElement)
    return tl.attrs[bn.ID]
end

function copy(x::TraitLikelihoodXMLElement)
    return TraitLikelihoodXMLElement(
        x.el,
        x.mbd_el,
        x.treeModel_el,
        x.extension_el,
        x.μ,
        x.pss,
        x.attrs,
        x.xml_name
    )
end

function make_xml(tl_el::TraitLikelihoodXMLElement)
    make_xml(tl_el.mbd_el)
    make_xml(tl_el.treeModel_el)
    make_xml(tl_el.extension_el)

    el = new_element(tl_el.xml_name)
    set_attributes(el, tl_el.attrs)

    add_ref_el(el, tl_el.mbd_el.el)
    add_ref_el(el, tl_el.treeModel_el.el)
    if !isnothing(tl_el.extension_el)
        add_ref_el(el, tl_el.extension_el)
    else
        tp_el = new_child(el, bn.TRAIT_PARAMETER)
        add_parameter_id(tp_el, bn.TRAIT_PARAMETER)
    end
    add_conjugate_root_prior(el, tl_el.μ, tl_el.pss)
    tl_el.el = el
    return el
end



function add_conjugate_root_prior(pel::XMLElement, μ::Vector{Float64}, pss::Float64)
    el = new_element(bn.CONJUGATE_ROOT_PRIOR)
    mean_el = new_element(bn.MEAN_PARAMETER)
    add_parameter(mean_el, value = μ)
    add_child(el, mean_el)

    pss_el = new_element(bn.PSS)
    add_parameter(pss_el, dim = 1, value = [pss])
    add_child(el, pss_el)

    add_child(pel, el)
    return el
end

function no_standardization!(el::TraitLikelihoodXMLElement)
    el.attrs[bn.SCALE_BY_TIME] = bn.FALSE
    el.attrs[bn.STANDARDIZE] = bn.FALSE
end

const DEFAULT_TRAITLIKELIHOOD_ATTRS = Dict(
        [(bn.ID, bn.TRAIT_DATA_LIKELIHOOD_ID),
        (bn.TRAIT_NAME, bn.DEFAULT_TRAIT_NAME),
        (bn.CACHE_BRANCHES, bn.TRUE),
        (bn.ALLOW_IDENTICAL, bn.TRUE),
        (bn.USE_TREE_LENGTH, bn.FALSE),
        (bn.SCALE_BY_TIME, bn.TRUE),
        (bn.REPORT_AS_MULTIVARIATE, bn.TRUE),
        (bn.INTEGRATE_INTERNAL_TRAITS, bn.TRUE),
        (bn.STANDARDIZE, bn.TRUE),
        (bn.ALLOW_SINGULAR, bn.FALSE)]
        )
