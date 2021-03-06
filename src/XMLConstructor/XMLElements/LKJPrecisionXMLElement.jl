mutable struct LKJPrecisionXMLElement <: MyXMLElement
    el::XMLOrNothing
    var::AbstractMatrix{Float64}
    diag_param::Parameter
    off_diag_parm::Parameter
    id::String

    function LKJPrecisionXMLElement(var_mat::AbstractArray{Float64},
                                    id::AbstractString)
        check_posdef(var_mat)
        vars, corrs = var_corrs(var_mat)
        diag_param = Parameter(vars, "$id.diagonalVar")
        diag_param.lower = 0.0
        off_diag_param = Parameter(corrs, "$id.correlation")
        off_diag_param.lower = -1.0
        off_diag_param.upper = 1.0
        return new(nothing, var_mat, diag_param, off_diag_param, id)
    end
end

function LKJPrecisionXMLElement(p::Int)
    return LKJPrecisionXMLElement(Diagonal(ones(p)), bn.DIFF_PREC_ID)
end


function make_xml(lkj::LKJPrecisionXMLElement)
    var, corrs = var_corrs(lkj.var)

    set_value(lkj.diag_param, var)
    set_value(lkj.off_diag_parm, corrs)

    el = new_element(bn.CACHED_MATRIX_INVERSE)
    set_id!(el, lkj.id)

    cs_el = new_child(el, bn.COMPOUND_SYMMETRIC_MATRIX)
    set_attributes(cs_el,
                   Dict(bn.AS_CORRELATION => bn.TRUE, bn.IS_CHOLESKY => bn.TRUE)
                  )

    diag_el = new_child(cs_el, bn.DIAGONAL)
    dp_el = make_xml(lkj.diag_param)
    add_child(diag_el, dp_el)
    od_el = new_child(cs_el, bn.OFF_DIAGONAL)
    op_el = make_xml(lkj.off_diag_parm)
    add_child(od_el, op_el)

    lkj.el = el
    return el
end

function var_corrs(var::AbstractArray{Float64})
    vars = diag(var)
    corr_mat = cov2corr(var)
    p = size(var, 1)
    n_corrs = div(p * (p - 1), 2)
    corrs = zeros(n_corrs)
    ind = 0
    for i = 1:p
        for j = (i + 1):p
            ind += 1
            corrs[ind] = corr_mat[i, j]
        end
    end
    return vars, corrs
end


################################################################################
## Priors
################################################################################

mutable struct LKJCorrelationPrior <: MyXMLElement
    el::XMLOrNothing
    corr_param::Parameter
    shape::Float64

    function LKJCorrelationPrior(corr_param::Parameter)
        return new(nothing, corr_param, 1.0)
    end
end

function make_xml(lkj::LKJCorrelationPrior)
    make_xml(lkj.corr_param)
    el = new_element(bn.LJK_CORRELATION_PRIOR)
    set_id!(el, get_id(lkj.corr_param.el) * ".prior")
    set_attribute(el, bn.SHAPE_PARAMETER, lkj.shape)
    dim = div(1 + isqrt(1 + 8 * length(lkj.corr_param)), 2)
    set_attribute(el, bn.DIMENSION, dim)
    set_attribute(el, bn.CHOLESKY, bn.TRUE)
    data_el = new_child(el, bn.DATA)
    add_ref_el(data_el, lkj.corr_param.el)
    lkj.el = el
    return el
end

mutable struct LogNormalPriorXMLElement <: MyXMLElement
    el::XMLOrNothing
    param::Parameter
    μ::Float64
    σ::Float64
    offset::Float64

    function LogNormalPriorXMLElement(param::Parameter)
        return new(nothing, param, 0.0, 1.0, 0.0)
    end
end

function make_xml(lnp::LogNormalPriorXMLElement)
    make_xml(lnp.param)
    el = new_element(bn.LOG_NORMAL_PRIOR)
    set_id!(el, get_id(lnp.param.el) * ".prior")
    set_attributes(el, Dict(
                            bn.MEAN => string(lnp.μ),
                            bn.STDEV => string(lnp.σ),
                            bn.OFFSET => string(lnp.offset),
                            bn.MEAN_IN_REAL_SPACE => bn.FALSE
                           )
                  )
    add_ref_el(el, lnp.param.el)
    lnp.el = el
    return el
end


function LKJPrecisionPriors(lkj::LKJPrecisionXMLElement)
    corr_prior = LKJCorrelationPrior(lkj.off_diag_parm)
    var_prior = LogNormalPriorXMLElement(lkj.diag_param)
    return CompoundXMLElement([corr_prior, var_prior])
end



