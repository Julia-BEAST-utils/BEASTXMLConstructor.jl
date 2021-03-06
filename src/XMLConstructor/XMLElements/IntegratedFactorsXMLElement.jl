abstract type LoadingsPriorXMLElement <: MyXMLElement end


mutable struct IntegratedFactorsXMLElement <: ModelExtensionXMLElement
    el::XMLOrNothing
    precision_prior_el::XMLOrNothing
    loadings::MyXMLElement
    loadings_prior::MyXMLElement
    precision::MyXMLElement
    precision_scale::Float64
    precision_shape::Float64
    treeModel::TreeModelXMLElement
    tree_trait_ind::Int
    standardize_traits::Bool
    id::String
end

function IntegratedFactorsXMLElement(treeModel_el::TreeModelXMLElement,
                                     k::Int; trait_ind::Int = 1,
                                     orthonormal::Bool = false)


    p = treeModel_el.trait_dims[trait_ind]

    L = default_loadings(k, p)
    if orthonormal
        load_param = ScaledOrthogonalMatrix(L, "L", "scale", "U")
    else
        load_param = MatrixParameter(L, "L", ["L$i" for i = 1:size(L, 1)])
    end

    prior = NormalDistributionLikelihood(load_param, "L.prior")


    precision = Parameter(ones(p), bn.FACTOR_PRECISION, lower=0.0)
    precision_scale = 1.0
    precision_shape = 1.0
    return IntegratedFactorsXMLElement(nothing, nothing,
                                load_param,
                                prior,
                                precision, precision_scale,
                                precision_shape,
                                treeModel_el,
                                trait_ind, true,
                                bn.DEFAULT_IF_NAME)
end

function name(::IntegratedFactorsXMLElement)
    return bn.INTEGRATED_FACTORS
end

function make_xmlelement(model::IntegratedFactorModel, tm::TreeModelXMLElement;
                         ind::Int = 1)
    ifx = IntegratedFactorsXMLElement(tm, tip_dimension(model), trait_ind = ind)
    set_loadings(ifx, model.L)
    set_precision(ifx, 1 ./ model.λ)
    return ifx
end

function copy(x::IntegratedFactorsXMLElement)
    return IntegratedFactorsXMLElement(
        x.el,
        x.precision_prior_el,
        x.loadings,
        x.loadings_prior,
        x.precision,
        x.precision_scale,
        x.precision_shape,
        x.treeModel,
        x.tree_trait_ind,
        x.standardize_traits,
        x.id
    )
end

function default_loadings(k::Int, p::Int)
    @assert k <= p
    L = zeros(k, p)
    for i = 1:k
        L[i, i] = 1.0
    end
    # L = randn(k, p)
    return L
end

function set_loadings!(ifxml::IntegratedFactorsXMLElement, L::S) where S <: AbstractArray{Float64, 2}
    set_mat!(ifxml.loadings, L)
end

function get_loadings_param(ifxml::IntegratedFactorsXMLElement)
    return ifxml.loadings
end

function get_precision(ifxml::IntegratedFactorsXMLElement)
    return ifxml.precision
end

function make_xml(ifxml::IntegratedFactorsXMLElement;
                    reference_precision::Bool = false)

    make_xml(ifxml.loadings)
    # ifxml.loadings_el = make_loadings(ifxml.loadings)
    actual_prior = nothing
    make_xml(ifxml.loadings_prior)

    el = new_element(bn.INTEGRATED_FACTORS)
    attrs = [(bn.ID, ifxml.id),
            (bn.TRAIT_NAME, ifxml.treeModel.node_traits[ifxml.tree_trait_ind])]
    set_attributes(el, attrs)

    if ifxml.standardize_traits
        set_attribute(el, bn.STANDARDIZE, bn.TRUE)
    else
        set_attribute(el, bn.STANDARDIZE, bn.FALSE)
    end

    l_el = new_child(el, bn.LOADINGS)
    add_ref_el(l_el, ifxml.loadings.el)

    p_el = new_child(el, bn.PRECISION)
    if reference_precision
        add_ref_el(p_el, ifxml.precision)
    else
        prec_el = make_xml(ifxml.precision)
        add_child(p_el, prec_el)
    end
    add_ref_el(el, ifxml.treeModel)
    tp_el = new_child(el, bn.TRAIT_PARAMETER)
    tp_id = ifxml.treeModel.param_names[ifxml.tree_trait_ind]
    add_ref_el(tp_el, bn.PARAMETER, tp_id)



    ifxml.el = el

    if !reference_precision
        ifxml.precision_prior_el = make_precision_prior(ifxml)
    end


    return el
end


function make_precision_prior(ifxml::IntegratedFactorsXMLElement)
    el = new_element(bn.GAMMA_PRIOR)
    set_attributes(el, [bn.ID => "$(bn.FACTOR_PRECISION).prior",
                        bn.SCALE => string(ifxml.precision_scale),
                        bn.SHAPE => string(ifxml.precision_shape)])
    add_ref_el(el, ifxml.precision)
    return el
end


function get_hmc_parameter(ifxml::IntegratedFactorsXMLElement)
    make_xml(ifxml)
    if typeof(ifxml.loadings) <: ScaledOrthogonalMatrix
        return ifxml.loadings.U
    end
    return ifxml.loadings.el
end

function make_loadings(L::AbstractArray{Float64, 2})
    el = new_element(bn.MATRIX_PARAMETER)
    set_attribute(el, bn.ID, bn.DEFAULT_LOADINGS_ID)
    n, p = size(L)
    for i = 1:n
        add_parameter(el, id = "$(bn.DEFAULT_LOADINGS_ID)$i",
                        value = L[i, :])
    end
    return el
end


function get_priors(xml::IntegratedFactorsXMLElement)
    make_xml(xml)

    return [xml.precision_prior_el; get_priors(xml.loadings_prior)]
end

function get_loggables(xml::IntegratedFactorsXMLElement)
    make_xml(xml)

    loggables = [get_loggables(xml.loadings); xml.precision]
    loggables = [loggables; get_loggables(xml.loadings_prior)]

    return loggables
end

function get_precision_prior(xml::IntegratedFactorsXMLElement)
    return xml.precision_prior_el
end

function set_shrinkage_mults!(ifxml::IntegratedFactorsXMLElement;
                        shapes::Vector{Float64} = Float64[],
                        scales::Vector{Float64} = Float64[],
                        set_scale::Bool = true,
                        set_mults::Bool = true)
    # msls = ifxml.msls
    # if isnothing(msls)
    #     error("No shrinkage prior on the integrated factors model.")
    # else
    #     if length(shapes) > 0
    #         fill_shrinkage_array!(msls.shapes, shapes)
    #     end
    #     if length(scales) > 0
    #         fill_shrinkage_array!(msls.scales, scales)
    #     end
    # end
    scale = set_shrinkage_mults!(ifxml.loadings_prior, shapes=shapes,
                                 scales=scales,
                                 set_mults = set_mults)
    if set_scale
        set_scale!(ifxml.loadings, scale)
    end
end

function fill_shrinkage_array!(to_fill::Vector{Float64},
                               fill_from::Vector{Float64})
    n = length(to_fill)
    m = length(fill_from)
    if n == m
        to_fill .= fill_from
    elseif (n - 1) == m
        to_fill[2:end] .= fill_from
    else
        error("Not implemented for these array dimensions.")
    end
end

function set_loadings(pfa::IntegratedFactorsXMLElement,
                      L::AbstractArray{Float64, 2})
    mat_param = pfa.loadings
    set_mat!(mat_param, L)
end


function set_precision(ifa::IntegratedFactorsXMLElement,
                       mat::AbstractArray{Float64, 2})
    if !isdiag(mat)
        error("Factor model must have diagonal precision")
    end
    set_precision(ifa, diag(mat))
end

function set_precision(ifa::IntegratedFactorsXMLElement,
                                precs::AbstractArray{Float64, 1})
    p = length(precs)
    if length(ifa.precision) != p
        error("Cannot set precision. Incompatible dimensions.")
    end
    for i = 1:p
        if precs[i] <= 0
            error("Precisions must be positive.")
        end
    end
    ifa.precision = precs
end

function get_loadings_scale(ifm::IntegratedFactorsXMLElement)
    loadings = ifm.loadings
    if typeof(loadings) <: ScaledOrthogonalMatrix
        return loadings.scale
    else
        error("Scale matrix does not exist " *
              "(loadings must be of type ScaledOrthogonalMatrix).")
    end
end

function get_loadings(ifm::IntegratedFactorsXMLElement)
    loadings = ifm.loadings
    if typeof(loadings) <: ScaledOrthogonalMatrix
        return loadings.U
    end
    return loadings
end
