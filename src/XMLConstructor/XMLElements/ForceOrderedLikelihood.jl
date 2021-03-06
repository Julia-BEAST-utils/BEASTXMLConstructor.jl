mutable struct ForceOrderedLikelihood <: MyXMLElement
    el::XMLOrNothing
    parameter::MyXMLElement

    function ForceOrderedLikelihood(p::MyXMLElement)
        return new(nothing, p)
    end
end

function make_xml(fol::ForceOrderedLikelihood)
    el = new_element("forceOrderedLikelihood")
    make_xml(fol.parameter)
    param_id = get_id(fol.parameter.el)
    add_ref_el(el, fol.parameter.el)
    set_id!(el, param_id * ".forceDescending")
    fol.el = el
    return el
end

function get_priors(fol::ForceOrderedLikelihood)
    make_xml(fol)
    return [fol.el]
end

mutable struct ColumnSwapOperator <: OperatorXMLElement
    el::XMLOrNothing
    matparam::MatrixParameter
    weight::Float64

    function ColumnSwapOperator(mp::MatrixParameter, weight::Float64)
        return new(nothing, mp, weight)
    end
end

function make_xml(cso::ColumnSwapOperator)
    el = new_element("columnSwapOperator")
    make_xml(cso.matparam)
    add_ref_el(el, cso.matparam.el)
    set_attribute(el, bn.WEIGHT, cso.weight)
    cso.el = el
    return el
end