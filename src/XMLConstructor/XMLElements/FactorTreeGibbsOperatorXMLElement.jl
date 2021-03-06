mutable struct FactorTreeGibbsOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    lfm::LatentFactorModelXMLElement
    tdl::TraitLikelihoodXMLElement
    weight::Real

    function FactorTreeGibbsOperatorXMLElement(lfm::LatentFactorModelXMLElement,
                    tdl::TraitLikelihoodXMLElement)
        return new(nothing, lfm, tdl, 1.0)
    end
end

function make_xml(op::FactorTreeGibbsOperatorXMLElement)
    el = new_element(bn.FACTOR_TREE_GIBBS_OP)
    set_attribute(el, bn.WEIGHT, op.weight)
    make_xml(op.lfm)
    make_xml(op.tdl)
    add_ref_el(el, op.lfm.el)
    add_ref_el(el, op.tdl.el)
    op.el = el
    return el
end
