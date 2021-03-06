mutable struct LoggablesXMLElement <: MyXMLElement
    els::Vector{<:MyXMLElement}
    already_made::Vector{Bool}
    function LoggablesXMLElement(els::Vector{<:MyXMLElement},
                                 already_made::Vector{Bool})
        if length(els) != length(already_made)
            error("Both arguments must be of the same length")
        end
        return new(els, already_made)
    end
end

function LoggablesXMLElement()
    return LoggablesXMLElement(MyXMLElement[], Bool[])
end

function LoggablesXMLElement(els::Vector{<:MyXMLElement};
                             already_made::Bool = true)
    return LoggablesXMLElement(els, fill(already_made, length(els)))
end

function make_xml(lg_el::LoggablesXMLElement)
    for el in lg_el.els
        make_xml(el)
    end
end

function join_loggables(lg1::LoggablesXMLElement, lg2::LoggablesXMLElement)
    return LoggablesXMLElement([lg1.els; lg2.els], [lg1.already_made; lg2.already_made])
end

import Base: push!
function push!(lg1, lg2)
    lg1.els = [lg1.els; lg2.els]
    lg1.already_made = [lg1.already_made; lg2.already_made]
end

function add_loggable(lg_el::LoggablesXMLElement, my_el::MyXMLElement;
                    already_made::Bool = false)
    push!(lg_el.els, my_el)
    push!(lg_el.already_made, already_made)
end


mutable struct MatrixInverseXMLElement <: SimpleXMLElement
    el::XMLOrNothing
    matrixParam_provider::MyXMLElement

    MatrixInverseXMLElement(my_xml::MyXMLElement) = new(nothing, my_xml)
end

function make_xml(mi_el::MatrixInverseXMLElement)
    make_xml(mi_el.matrixParam_provider)
    mp_el = get_matrix_parameter(mi_el.matrixParam_provider)

    el = new_element(bn.MATRIX_INVERSE)
    id = attribute(mp_el, bn.ID)
    set_attribute(el, bn.ID, "inverse.$id")
    add_ref_el(el, mp_el)

    mi_el.el = el
    return el
end

function get_loggables(mi_el::MatrixInverseXMLElement)
    make_xml(mi_el)
    return [mi_el.el]
end


mutable struct CorrelationMatrixXMLElement <: SimpleXMLElement
    el::XMLOrNothing
    matrixParam_provider::MyXMLElement
    invert::Bool

    CorrelationMatrixXMLElement(my_xml::MyXMLElement) = new(nothing, my_xml, false)
    CorrelationMatrixXMLElement(my_xml::MyXMLElement, invert::Bool) = new(nothing, my_xml, invert)
end

function make_xml(cm_el::CorrelationMatrixXMLElement)
    make_xml(cm_el.matrixParam_provider)
    mp_el = get_matrix_parameter(cm_el.matrixParam_provider)

    el = new_element(bn.CORRELATION_MATRIX)
    id = attribute(mp_el, bn.ID)
    nm = "correlation"
    if cm_el.invert
        nm = "$nm.inverted"
    end
    set_attribute(el, bn.ID, "$nm.$id")
    set_attribute(el, bn.INVERT, cm_el.invert)
    add_ref_el(el, mp_el)

    cm_el.el = el
    return el
end

function get_loggables(cm_el::CorrelationMatrixXMLElement)
    make_xml(cm_el)
    return [cm_el.el]
end

mutable struct VarianceProportionXMLElement <: SimpleXMLElement
    el::XMLOrNothing
    td_el::TraitLikelihoodXMLElement
    tm_el::TreeModelXMLElement
    rm_el::RepeatedMeasuresXMLElement
    mbd_el::MBDXMLElement
    id::String
    emp_var::Bool

    function VarianceProportionXMLElement(td::TraitLikelihoodXMLElement,
                                tm::TreeModelXMLElement,
                                rm::RepeatedMeasuresXMLElement,
                                mbd::MBDXMLElement)

        return new(nothing, td, tm, rm, mbd, bn.VARIANCE_PROPORTION_STATISTIC,
                    false)
    end

end

function make_xml(vp_el::VarianceProportionXMLElement)
    make_xml(vp_el.td_el)
    make_xml(vp_el.tm_el)
    make_xml(vp_el.rm_el)
    make_xml(vp_el.mbd_el)

    el = new_element(bn.VARIANCE_PROPORTION_STATISTIC)
    set_attribute(el, bn.ID, vp_el.id)
    set_attribute(el, bn.MATRIX_RATIO, bn.COHERITABILITY)
    set_attribute(el, bn.USE_EMPIRICAL_VARIANCE, vp_el.emp_var)

    add_ref_el(el, vp_el.td_el.el)
    add_ref_el(el, vp_el.tm_el.el)
    add_ref_el(el, vp_el.rm_el.el)
    add_ref_el(el, vp_el.mbd_el.el)

    vp_el.el = el
    return el
end

function get_loggables(vp_el::VarianceProportionXMLElement)
    make_xml(vp_el)
    return [vp_el.el]
end

mutable struct ModelExtensionLoggerXMLElement <: MyXMLElement
    el::XMLOrNothing
    tl_el::TraitLikelihoodXMLElement
    ext_el::MyXMLElement

    ModelExtensionLoggerXMLElement(tl::TraitLikelihoodXMLElement, my::MyXMLElement) =
        new(nothing, tl, my)
end

function make_xml(mel_el::ModelExtensionLoggerXMLElement)
    el = new_element(bn.MODEL_EXTENSION_TRAIT_LOGGER)
    set_attributes(el, [(bn.ID, bn.MODEL_EXTENSION_TRAIT_LOGGER),
                        (bn.DIMENSIONS, bn.MISSING),
                        (bn.TRAIT_NAME, mel_el.tl_el.attrs[bn.TRAIT_NAME])])

    make_xml(mel_el.tl_el)
    make_xml(mel_el.ext_el)

    add_ref_el(el, mel_el.tl_el.el)
    add_ref_el(el, mel_el.ext_el.el)

    mel_el.el = el
    return el
end

function get_loggables(mel_el::ModelExtensionLoggerXMLElement)
    return [mel_el.el]
end

# mutable struct MatrixValidationXMLElement <: MyXMLElement
#     el::XMLOrNothing
#     matrixParam_provider::MyXMLElement
#     true_parameter::AbstractArray{Float64}
#
#     MatrixValidationXMLElement(my_xml::MyXMLElement, x::AbstractArray{Float64, 2})
#         = new(my_xml, x)
# end
#
# function make_xml(mv_el::MatrixValidationXMLElement)
#     make_xml(mv_el.matrixParam_provider)
#
#     el = new_element(bn.MATRIX_VALIDATION)
#     tp_el = new_child(el, bn.TRUE_PARAMETER)
#     add_matrix_parameter(tp_el, true_parameter)
#     ip_el = new_child(el, bn.INFERRED_PARAMETER)
#     mat_param = get_matrix_parameter(mv_el.matrixParam_provider)
#     add_ref_el(ip_el, mat_param)
#
