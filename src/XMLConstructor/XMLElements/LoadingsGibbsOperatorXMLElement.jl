mutable struct LoadingsGibbsOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    if_xml::IntegratedFactorsXMLElement
    tdl_xml::TraitLikelihoodXMLElement
    sparsity::String
    weight::Real

    function LoadingsGibbsOperatorXMLElement(
            if_xml::IntegratedFactorsXMLElement,
            tdl_xml::TraitLikelihoodXMLElement)

        return new(nothing, if_xml, tdl_xml, bn.NONE, 1.0)
    end
end

function make_xml(lgoxml::LoadingsGibbsOperatorXMLElement)
    el = new_element(bn.LOADINGS_GIBBS_OP)
    set_attribute(el, bn.WEIGHT, lgoxml.weight)
    set_attributes(el, [bn.RANDOM_SCAN => bn.FALSE,
                        bn.NEW_MODE => bn.TRUE,
                        bn.CONSTRAINT => bn.NONE,
                        bn.SPARSITY => lgoxml.sparsity])

    make_xml(lgoxml.if_xml)
    make_xml(lgoxml.tdl_xml)

    add_ref_el(el, lgoxml.if_xml.el)
    add_ref_el(el, lgoxml.tdl_xml.el)
    add_ref_el(el, get_normal_prior(lgoxml.if_xml.loadings_prior))

    lgoxml.el = el
end

function sparsity_constraint!(lgo::LoadingsGibbsOperatorXMLElement, constraint::String)
    lgo.sparsity = constraint
end

function get_parameter(op::LoadingsGibbsOperatorXMLElement)
    return get_loadings(op.if_xml)
end