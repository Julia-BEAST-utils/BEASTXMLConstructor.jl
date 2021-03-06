mutable struct MBDXMLElement <: MyXMLElement
    el::XMLOrNothing
    precision::MyXMLElement
    precision_prior::MyXMLElement
    is_random::Bool
end

function MBDXMLElement(precision::AbstractArray{Float64, 2},
                       prior_scale::AbstractArray{Float64, 2},
                       is_random::Bool,
                       diagonal_prec::Bool)
    if diagonal_prec
        if !isdiag(prior_scale)
            error("Matrix must be diagonal.")
        else
            precision_param = DiagonalMatrixParameter(diag(prior_scale))
            prior = NothingXMLElement()
        end
    else
        precision_param = MatrixParameter(precision, bn.DIFF_PREC_ID)
        prior = WishartPriorXMLElement(prior_scale, precision_param)
    end
    if !is_random && diagonal_prec
        error("Not implemented") # TODO: gamma prior?
    end
    return MBDXMLElement(nothing, precision_param, prior, is_random)
end


function MBDXMLElement(p::Int; diagonal_prec::Bool = false)
    return MBDXMLElement(Diagonal(ones(p)),
                         Diagonal(ones(p)),
                         true,
                         diagonal_prec)
end

function MBDXMLElement(dm::DiffusionModel)
    P = inv(dm.Σ)
    return MBDXMLElement(P, Diagonal(ones(size(P, 1))), true, false)
end


function make_xml(ml::MBDXMLElement)
    el = new_element(bn.MULTIVARIATE_DIFFUSION_MODEL)
    set_attribute(el, bn.ID, bn.DIFFUSION_ID)
    prec_el = new_element(bn.PRECISION_MATRIX)

    mat_el = make_xml(ml.precision)
    add_child(prec_el, mat_el)

    add_child(el, prec_el)
    ml.el = el

    if ml.is_random
        make_xml(ml.precision_prior)
    end
end



function get_loggables(el::MBDXMLElement)
    if !el.is_random
        return Vector{XMLElement}(undef, 0)
    end
    make_xml(el)
    return [el.el[bn.PRECISION_MATRIX][1][bn.MATRIX_PARAMETER][1]]

end

function get_matrix_parameter(el::MBDXMLElement)
    return get_loggables(el)[1]
end

function get_priors(xml::MBDXMLElement)
    if !xml.is_random
        return Vector{XMLElement}(undef, 0)
    end
    make_xml(xml)
    return [xml.precision_prior.el]
end

function set_precision(mbd::MBDXMLElement, prec::AbstractArray{Float64, 2})
    check_posdef(prec)
    mbd.precision = prec
end


################################################################################
## Priors
################################################################################

