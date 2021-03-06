mutable struct WishartPriorXMLElement <: MyXMLElement
    el::XMLOrNothing
    df::Int
    scale_mat::AbstractMatrix{Float64}
    param::MatrixParameter
    id::String

    function WishartPriorXMLElement(scale_mat::AbstractMatrix{Float64}, df::Int,
                                    ref_param::MatrixParameter)
        check_posdef(scale_mat)
        p = size(scale_mat, 1)
        if df < p
            error("Degrees of freedom must be at least the matrix dimension.")
        end
        return new(nothing, df, scale_mat, ref_param, ref_param.id * ".prior")
    end
end

function WishartPriorXMLElement(param::MatrixParameter)
    p = size(param, 1)
    return WishartPriorXMLElement(Diagonal(ones(p)), p, param)
end

function WishartPriorXMLElement(mat::AbstractMatrix{Float64},
                                param::MatrixParameter)
    if size(mat) != size(param)
        error("Incompatible dimensions.")
    end
    return WishartPriorXMLElement(mat, size(mat, 1), param)
end

function make_xml(wp::WishartPriorXMLElement)
    make_xml(wp.param)
    el = new_element(bn.MULTIVARIATE_WISHART_PRIOR)
    set_attributes(el, [(bn.ID, wp.id), (bn.DF, string(wp.df))])
    scale_el = new_child(el, bn.SCALE_MATRIX)
    add_matrix_parameter(scale_el, wp.scale_mat)
    data_el = new_child(el, bn.DATA)
    add_ref_el(data_el, wp.param.el)
    wp.el = el
    return el
end



# function make_Wishart_prior(scale::AbstractArray{Float64, 2},
#                             mp_el::XMLElement,
#                             id::String)

#     el = new_element(bn.MULTIVARIATE_WISHART_PRIOR)
#     set_attributes(el, [(bn.ID, id), (bn.DF, string(size(scale, 1)))])
#     scale_el = new_child(el, bn.SCALE_MATRIX)
#     add_matrix_parameter(scale_el, scale)
#     data_el = new_child(el, bn.DATA)
#     add_ref_el(data_el, mp_el)
#     return el
# end