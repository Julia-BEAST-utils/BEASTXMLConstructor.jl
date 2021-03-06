mutable struct OldLoadingsGibbsOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    lfm::LatentFactorModelXMLElement
    weight::Real
    randomScan::Bool

    function OldLoadingsGibbsOperatorXMLElement(lfm::LatentFactorModelXMLElement)
        return new(nothing, lfm, 1.0, false)
    end

end

function make_xml(lgo::OldLoadingsGibbsOperatorXMLElement)
    el = new_element(bn.LOADINGS_GIBBS_OP)
    set_attribute(el, bn.WEIGHT, lgo.weight)
    set_attribute(el, bn.RANDOM_SCAN, lgo.randomScan)
    set_attribute(el, bn.NEW_MODE, bn.TRUE)
    set_attributes(el, [bn.SPARSITY => bn.NONE, bn.CONSTRAINT => bn.NONE])

    make_xml(lgo.lfm)
    add_ref_el(el, lgo.lfm.el)
    add_ref_el(el, make_xml(get_loadings_prior(lgo.lfm.loadings_prior)))

    lgo.el = el
    return el

end
