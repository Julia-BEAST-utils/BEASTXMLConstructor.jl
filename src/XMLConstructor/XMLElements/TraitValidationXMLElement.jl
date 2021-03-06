mutable struct TraitValidationXMLElement <: MyXMLElement
    el::XMLOrNothing
    treeModel_el::TreeModelXMLElement
    traitLikelihood_el::TraitLikelihoodXMLElement
    mis_ind::Int
    obs_ind::Int
    standardize::Bool
    id::String

    function TraitValidationXMLElement(tm_el::TreeModelXMLElement,
                    tl_el::TraitLikelihoodXMLElement)

        @assert length(tm_el.node_traits) == 2

        mis_ind = 1
        obs_ind = 2

        return new(nothing, tm_el, tl_el, mis_ind, obs_ind, false, bn.TRAIT_VALIDATION)
    end
end

mutable struct FactorValidationXMLElement <: MyXMLElement
    el::XMLOrNothing
    ifl::IntegratedFactorsXMLElement
    tdl::TraitLikelihoodXMLElement
    tml::TreeModelXMLElement
    trait_name::String
    trait_param_name::String
    id::String
end

function FactorValidationXMLElement(ifl::IntegratedFactorsXMLElement,
                                    tdl::TraitLikelihoodXMLElement,
                                    tml::TreeModelXMLElement
                                    )

    @assert length(tml.param_names) == 2

    return FactorValidationXMLElement(nothing, ifl, tdl, tml,
                                      tml.node_traits[2], tml.param_names[2],
                                      "factorValidation")
end

mutable struct CrossValidationXMLElement <: MyXMLElement
    el::XMLOrNothing
    validationProvider_el::MyXMLElement
    log_sum::Bool
    validation_type::String
    id::String

    CrossValidationXMLElement(xml::MyXMLElement) =
            new(nothing, xml, false, bn.SQUARED_ERROR, bn.CROSS_VALIDATION)
end

function set_validation_type!(cv::CrossValidationXMLElement, t::String)
    cv.validation_type = t
end

function set_log_sum!(cv::CrossValidationXMLElement, b::Bool)
    cv.log_sum = b
end


function make_xml(cv_el::CrossValidationXMLElement)

    el = new_element(bn.CROSS_VALIDATION)
    set_attributes(el, [(bn.ID, cv_el.id),
                        (bn.LOG_SUM, string(cv_el.log_sum)),
                        (bn.TYPE, cv_el.validation_type)])

    make_xml(cv_el.validationProvider_el)
    add_child(el, cv_el.validationProvider_el.el)

    cv_el.el = el
    return el
end


function make_xml(tv_el::TraitValidationXMLElement)

    obs_ind = tv_el.obs_ind
    mis_ind = tv_el.mis_ind

    el = new_element(bn.TRAIT_VALIDATION_PROVIDER)
    attrs = [(bn.ID, bn.TRAIT_VALIDATION),
            (bn.TRAIT_NAME, tv_el.treeModel_el.node_traits[obs_ind]),
            (bn.INFERRED_TRAIT, tv_el.traitLikelihood_el.attrs[bn.TRAIT_NAME]),
            (bn.STANDARDIZE, string(tv_el.standardize))]
    set_attributes(el, attrs)

    make_xml(tv_el.traitLikelihood_el)
    add_ref_el(el, tv_el.traitLikelihood_el.el)

    tp_el = new_child(el, bn.TRAIT_PARAMETER)
    add_ref_el(tp_el, bn.PARAMETER, tv_el.treeModel_el.param_names[obs_ind])

    tv_el.el = el
    return el
end

function make_xml(fac_val::FactorValidationXMLElement)
    el = new_element(bn.FACTOR_VALIDATION)

    make_xml(fac_val.ifl)
    make_xml(fac_val.tdl)
    make_xml(fac_val.tml)

    add_ref_el(el, fac_val.ifl.el)
    add_ref_el(el, fac_val.tdl.el)
    trait_name = attribute(fac_val.tdl.el, bn.TRAIT_NAME)
    set_attribute(el, bn.ID, fac_val.id)
    set_attribute(el, bn.TRAIT_NAME, trait_name)

    trait_provider = new_child(el, bn.TRAIT_DATA_PROVIDER)
    add_ref_el(trait_provider, fac_val.tml.el)
    set_attribute(trait_provider, bn.TRAIT_NAME, fac_val.trait_name)
    tp_el = new_child(trait_provider, bn.TRAIT_PARAMETER)
    add_ref_el(tp_el, bn.PARAMETER, fac_val.trait_param_name)

    fac_val.el = el
end


function get_loggables(cv_el::CrossValidationXMLElement)
    make_xml(cv_el)
    return [cv_el.el]
end
