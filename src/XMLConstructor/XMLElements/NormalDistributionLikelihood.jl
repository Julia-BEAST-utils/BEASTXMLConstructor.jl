mutable struct NormalDistributionLikelihood <: MyXMLElement
    el::XMLOrNothing
    param::MyXMLElement
    mean::Float64
    sd::Float64
    id::String

    function NormalDistributionLikelihood(param::MyXMLElement, μ::Float64,
                                          σ::Float64, id::String)
        return new(nothing, param, μ, σ, id)
    end
end

function NormalDistributionLikelihood(param::MyXMLElement, id::String)
    return NormalDistributionLikelihood(param, 0.0, 1.0, id)
end



function make_xml(ndl::NormalDistributionLikelihood)
    el = new_element(bn.DISTRIBUTION_LIKELIHOOD)
    set_attribute(el, bn.ID, ndl.id)

    d_el = new_child(el, bn.DATA)
    make_xml(ndl.param)
    add_ref_el(d_el, ndl.param.el)
    dist_el = new_child(el, bn.DISTRIBUTION)

    ndm_el = new_child(dist_el, bn.NORMAL_DISTRIBUTION_MODEL)

    m_el = new_child(ndm_el, bn.MEAN)
    add_parameter(m_el, value = [ndl.mean])
    std_el = new_child(ndm_el, bn.STDEV)
    add_parameter(std_el, value =[ndl.sd], lower=0.0)

    ndl.el = el

    return el
end

function get_normal_prior(ndl::NormalDistributionLikelihood)
    make_xml(ndl)
    return ndl.el
end

function get_priors(ndl::NormalDistributionLikelihood)
    make_xml(ndl)
    return ndl.el
end

function get_loggables(::NormalDistributionLikelihood)
    return MyXMLElement[]
end

function make_loadings_gradient(loadings::MyXMLElement,
                                ndl::NormalDistributionLikelihood)
    el = new_element(bn.GRADIENT)
    make_xml(ndl)
    make_xml(loadings)

    add_ref_el(el, ndl.el)
    add_ref_el(el, loadings.el)

    return el
end

function get_loadings_prior(ndl::NormalDistributionLikelihood)
    return ndl
end

