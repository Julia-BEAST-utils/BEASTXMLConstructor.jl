mutable struct IntegratedFactorsGibbsOperator <: OperatorXMLElement
    el::XMLOrNothing
    sampled_tree_likelihood::TraitLikelihoodXMLElement
    integrated_tree_likelihood::TraitLikelihoodXMLElement
    factor_likelihood::IntegratedFactorsXMLElement
    weight::Float64

    function IntegratedFactorsGibbsOperator(
                    sampled_tree_likelihood::TraitLikelihoodXMLElement,
                    tree_likelihood::TraitLikelihoodXMLElement,
                    factor_likelihood::IntegratedFactorsXMLElement,
                    weight::Float64)
        return new(nothing, sampled_tree_likelihood,
                    tree_likelihood, factor_likelihood, weight)
    end
end

function make_xml(ifg::IntegratedFactorsGibbsOperator)
    el = new_element(bn.INTEGRATED_FACTOR_GIBBS)
    set_attribute(el, bn.WEIGHT, ifg.weight)

    sampled_el = make_xml(ifg.sampled_tree_likelihood)
    tp_el = find_element(sampled_el, bn.TRAIT_PARAMETER)
    fac_mat_el = find_element(tp_el, bn.PARAMETER)
    add_ref_el(el, fac_mat_el)

    integrated_el = make_xml(ifg.integrated_tree_likelihood)
    add_ref_el(el, integrated_el)

    fac_el = make_xml(ifg.factor_likelihood)
    add_ref_el(el, fac_el)

    ifg.el = el
    return el
end
