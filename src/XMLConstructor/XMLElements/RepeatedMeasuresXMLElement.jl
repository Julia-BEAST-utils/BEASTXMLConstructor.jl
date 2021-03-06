mutable struct RepeatedMeasuresXMLElement <: ModelExtensionXMLElement
    el::XMLOrNothing
    precision::MatrixParameter
    precision_prior::WishartPriorXMLElement
    treeModel_el::TreeModelXMLElement
    tree_trait_ind::Int
    trait_name::String
    tp_name::String
    standardize_traits::Bool


    function RepeatedMeasuresXMLElement(treeModel_el::TreeModelXMLElement;
                                        trait_ind::Int = 1)

        p = treeModel_el.trait_dims[trait_ind]
        trait_name = treeModel_el.node_traits[trait_ind]
        tp_name = treeModel_el.param_names[trait_ind]
        precision = MatrixParameter(Diagonal(ones(p)), bn.DEFAULT_RM_PREC_NAME)
        wishart_prior = WishartPriorXMLElement(precision)
        return new(nothing, precision, wishart_prior,
                   treeModel_el,
                   trait_ind, trait_name, tp_name, true)
    end
end

function make_xmlelement(rm::ResidualVarianceModel, tm::TreeModelXMLElement;
                         ind::Int = 1)
    rmxml = RepeatedMeasuresXMLElement(tm, trait_ind = ind)
    set_precision(rmxml, inv(rm.Γ))
    return rmxml
end

function make_xml(rm_el::RepeatedMeasuresXMLElement)

    # repeatedMeasures xml element
    el = new_element(bn.REPEATED_MEASURES)
    attrs = [(bn.ID, bn.DEFAULT_RM_NAME), (bn.TRAIT_NAME, rm_el.trait_name)]
    if rm_el.standardize_traits
        set_attribute(el, bn.STANDARDIZE, bn.TRUE)
    end
    set_attributes(el, attrs)
    treeModel_el = make_xml(rm_el.treeModel_el)
    add_ref_el(el, treeModel_el)
    tp_el = new_child(el, bn.TRAIT_PARAMETER)
    add_ref_el(tp_el, bn.PARAMETER, rm_el.tp_name)
    sp_el = new_child(el, bn.SAMPLING_PRECISION)
    mp_el = make_xml(rm_el.precision)
    add_child(sp_el, mp_el)

    # make_xml(rm_el.mbd_el)

    # add_ref_el(el, rm_el.mbd_el.el)
    rm_el.el = el

    # prior xml element

    make_xml(rm_el.precision_prior)

    return el
end

function get_loggables(el::RepeatedMeasuresXMLElement)
    make_xml(el)
    return [el.el[bn.SAMPLING_PRECISION][1][bn.MATRIX_PARAMETER][1]]
end

function get_matrix_parameter(el::RepeatedMeasuresXMLElement)
    return get_loggables(el)[1]
end

function get_priors(xml::RepeatedMeasuresXMLElement)
    make_xml(xml)
    return [xml.precision_prior.el]
end

function get_precision_prior(xml::RepeatedMeasuresXMLElement)
    return xml.precision_prior.el
end

function set_precision(rm::RepeatedMeasuresXMLElement,
                       mat::AbstractArray{Float64, 2})
    check_posdef(mat)
    set_mat!(rm.precision, mat)
end

