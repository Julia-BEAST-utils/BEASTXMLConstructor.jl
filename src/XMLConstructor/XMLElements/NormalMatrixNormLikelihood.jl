mutable struct NormalMatrixNormLikelihood <: MyXMLElement
    el::XMLOrNothing
    precision_parameter::MyXMLElement
    matrix_parameter::MyXMLElement
    id::String

    function NormalMatrixNormLikelihood(prec::MyXMLElement, mat::MyXMLElement,
                                        id::String)
        return new(nothing, prec, mat, id)
    end
end


function make_xml(nmnl::NormalMatrixNormLikelihood)
    @unpack el, precision_parameter, matrix_parameter, id = nmnl
    el = new_element(bn.NORMAL_MATRIX_NORM)
    set_id!(el, id)

    gp_el = new_child(el, bn.GLOBAL_PRECISION)
    make_xml(precision_parameter)
    add_ref_el(gp_el, precision_parameter.el)

    m_el = new_child(el, bn.MATRIX)
    make_xml(matrix_parameter)
    add_ref_el(m_el, matrix_parameter.el)

    nmnl.el = el
    return el
end

function get_normal_prior(nmnl::NormalMatrixNormLikelihood)
    make_xml(nmnl)
    return nmnl.el
end
