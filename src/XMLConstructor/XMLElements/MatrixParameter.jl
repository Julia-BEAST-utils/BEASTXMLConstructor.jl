function add_ref_el(el::XMLElement, param::MyXMLElement)
    ref_el = new_element(name(param))
    set_attribute(ref_el, bn.IDREF, get_id(param))
    add_child(el, ref_el)
end




################################################################################
## Parameter
################################################################################
import Base: size

mutable struct Parameter <: MyXMLElement
    el::XMLOrNothing
    val::T where T <: AbstractVector{Float64}
    dim::Int
    lower::Float64
    upper::Float64
    id::String
end

function Parameter(val::Vector{Float64}, id::String; lower::Float64 = NaN)
    return Parameter(nothing, val, length(val), lower, NaN, id)
end

function Parameter(val::Vector{<:Real})
    return Parameter(val, "")
end


function make_parameter(;id::String = "",
        value::AbstractArray{T} = Vector{Int8}(undef, 0),
        lower::Float64 = NaN, upper::Float64 = NaN,
        dim::Int = length(value)) where T <: Real

    el = new_element(bn.PARAMETER)
    if id != ""
        set_attribute(el, bn.ID, id)
    end
    if length(value) != 0
        set_attribute(el, bn.VALUE, join(value, ' '))
    end
    if !isnan(lower)
        set_attribute(el, bn.LOWER, lower)
    end
    if !isnan(upper)
        set_attribute(el, bn.UPPER, upper)
    end
    if dim != length(value)
        set_attribute(el, bn.DIMENSION, dim)
    end

    return el
end

function name(::Parameter)
    return bn.PARAMETER
end

function add_parameter(pel::XMLElement; id::String = "",
                       value::AbstractArray{<:Real} = Float64[],
                       lower::Float64 = NaN, upper::Float64 = NaN,
                       dim::Int = length(value)) where T <: Real

    el = make_parameter(id = id, value = value, lower = lower, upper = upper,
                        dim = dim)
    add_child(pel, el)
    return el
end


function make_xml(param::Parameter)
    el = make_parameter(id = param.id, value = param.val, lower = param.lower,
                        upper = param.upper, dim = param.dim)
    param.el = el
    return el
end

function set_value(p::Parameter, val::AbstractVector{Float64})
    if length(val) != length(p.val)
        error("Incompatible dimensions.")
    end
    p.val .= val
end

import Base: length
function length(p::Parameter)
    return p.dim
end

function get_loggables(p::Parameter)
    make_xml(p)
    return [p.el]
end


function get_hmc_parameter(param::Parameter)
    make_xml(param)
    return param.el
end


################################################################################
## MatrixParameter
################################################################################


mutable struct MatrixParameter <: MyXMLElement
    el::XMLOrNothing
    mat::AbstractMatrix{Float64}
    id::String
    ids::Vector{String}
end

function MatrixParameter(mat::AbstractMatrix{Float64}, id::String, ids::Vector{String})
    m, n = size(mat)
    if length(ids) != m
        error("The matrix has $m rows, but $(length(ids)) row ids were provided.")
    end

    return MatrixParameter(nothing, mat, id, ids)
end

function MatrixParameter(mat::AbstractMatrix{Float64}, id::String)
    m = size(mat, 1)
    return MatrixParameter(mat, id, ["" for i = 1:m])
end

function name(::MatrixParameter)
    return bn.MATRIX_PARAMETER
end

function make_xml(mp::MatrixParameter)
    el = new_element(bn.MATRIX_PARAMETER)
    set_attribute(el, bn.ID, mp.id)
    n = size(mp.mat, 1)
    for i = 1:n
        add_parameter(el, id = mp.ids[i], value = @view mp.mat[i, :])
    end
    mp.el = el
    return el
end

function set_mat!(mp::MatrixParameter, mat::S) where S <: AbstractMatrix{Float64}
    if size(mat) != size(mp.mat)
        error("Supplied matrix has dimensions $(size(mat)), while matrix " *
              "parameter has dimensions $(size(mp.mat))")
    end

    mp.mat = mat
end

function size(p::MatrixParameter)
    return size(p.mat)
end

function size(p::MatrixParameter, dim::Int)
    return size(p.mat, dim)
end

function get_loggables(mp::MatrixParameter)
    make_xml(mp)
    return [mp.el]
end

################################################################################
## DiagonalParameter
################################################################################

mutable struct DiagonalMatrixParameter <: MyXMLElement
    el::XMLOrNothing
    vals::AbstractArray{Float64}
    lower::Float64
    upper::Float64
    id::String
end

function DiagonalMatrixParameter(vals::AbstractVector{Float64}, id::String)
    return DiagonalMatrixParameter(nothing, vals, NaN, NaN, id)
end

function DiagonalMatrixParameter(vals::AbstractVector{Float64})
    return DiagonalMatrixParameter(nothing, vals, NaN, NaN, "")
end

function make_xml(dmp::DiagonalMatrixParameter)
    el = new_element(bn.DIAGONAL_MATRIX)
    if dmp.id != ""
        set_id!(el, dmp.id)
    end
    p_el = make_parameter(value = dmp.vals, lower=dmp.lower, upper=dmp.upper)
    add_child(el, p_el)
    dmp.el = el
    return el
end

function size(d::DiagonalMatrixParameter)
    p = length(d.vals)
    return (p, p)
end

function size(d::DiagonalMatrixParameter, dim::Int)
    p = length(d.vals)
    return p
end

################################################################################
## Multiplicative parameter
################################################################################

mutable struct MultiplicativeParameter <: MyXMLElement
    el::XMLOrNothing
    param::Parameter
    id::String

    function MultiplicativeParameter(param::Parameter, id::String)
        return new(nothing, param, id)
    end
end

function make_xml(mp::MultiplicativeParameter)
    el = new_element(bn.MULTIPLICATIVE_PARAMETER)
    set_id!(el, mp.id)

    make_xml(mp.param)
    add_ref_el(el, mp.param.el)

    mp.el = el
    return el
end

function name(::MultiplicativeParameter)
    return bn.MULTIPLICATIVE_PARAMETER
end

################################################################################
## Masked parameter
################################################################################

mutable struct MaskedParameter <: MyXMLElement
    el::XMLOrNothing
    param::MyXMLElement
    mask::Vector{Float64}

    function MaskedParameter(param::MyXMLElement, mask::Vector{Float64})
        return new(nothing, param, mask)
    end
end

function make_xml(mp::MaskedParameter)
    el = new_element(bn.MASKED_PARAMETER)
    add_ref_el(el, mp.param)
    m_el = new_child(el, bn.MASK)
    add_child(m_el, make_xml(Parameter(mp.mask)))
    mp.el = el
    return el
end


