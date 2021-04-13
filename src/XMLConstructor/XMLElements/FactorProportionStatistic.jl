mutable struct FactorProportionStatistic <: MyXMLElement
    el::XMLOrNothing
    factor_model::IntegratedFactorsXMLElement
    trait_model::TraitLikelihoodXMLElement
    id::String

    function FactorProportionStatistic(factor_model::IntegratedFactorsXMLElement,
                                        trait_model::TraitLikelihoodXMLElement,
                                        id::String)
        return new(nothing, factor_model, trait_model, id)
    end
end

function name(::FactorProportionStatistic)
    return bn.FACTOR_PROPORTION_STATISTIC
end

function make_xml(fps::FactorProportionStatistic)
    el = new_element(name(fps))

    set_attribute(el, bn.ID, fps.id)
    add_ref_el(el, fps.factor_model)
    add_ref_el(el, fps.trait_model)

    fps.el = el
    return el
end

function get_loggables(fps::FactorProportionStatistic)
    return fps
end