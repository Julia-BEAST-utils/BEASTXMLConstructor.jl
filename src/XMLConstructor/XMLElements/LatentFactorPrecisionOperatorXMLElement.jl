mutable struct LatentFactorModelPrecisionOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    lfm::LatentFactorModelXMLElement
    weight::Real

    function LatentFactorModelPrecisionOperatorXMLElement(
                        lfm::LatentFactorModelXMLElement)
        return new(nothing, lfm, 1.0)
    end
end

function make_xml(op::LatentFactorModelPrecisionOperatorXMLElement)
    el = new_element(bn.LAT_FAC_PREC_OP)
    set_attribute(el, bn.WEIGHT, op.weight)
    set_attribute(el, bn.RANDOM_SCAN, bn.FALSE)

    make_xml(op.lfm)
    add_ref_el(el, op.lfm.el)
    add_ref_el(el, op.lfm.precision_prior_el)
    op.el = el
    return el
end
