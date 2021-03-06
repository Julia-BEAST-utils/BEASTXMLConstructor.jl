abstract type Preconditioner end

mutable struct HMCOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    param_provider::MyXMLElement
    grads::Array{<:MyXMLElement}
    already_made::Vector{Bool}
    weight::Real
    n_steps::Int
    step_size::Float64
    draw_variance::Float64
    auto_optimize::Bool
    check_grad::Int
    grad_tolerance::Float64
    geodesic::Bool
    transform::String
    mask::Vector{<:Real}
    preconditioning::Preconditioner

    function HMCOperatorXMLElement(param_provider::MyXMLElement,
                        grads::Array{<:MyXMLElement};
                        geodesic::Bool = false,
                        already_made::Vector{Bool} = fill(false, length(grads)),
                        transform::String = "")
        return new(nothing, param_provider, grads, already_made,
                    1.0, 5, 0.05, 1.0, true, 0, 1e-3, geodesic, transform,
                    Float64[],
                    NoPreconditioning())
    end

end

function make_xml(hmcxml::HMCOperatorXMLElement)
    name = hmcxml.geodesic ? bn.GEODESIC_HMC : bn.HMC
    el = new_element(name)
    set_attributes(el, [bn.WEIGHT => string(hmcxml.weight),
                        bn.NSTEPS => string(hmcxml.n_steps),
                        bn.STEP_SIZE => string(hmcxml.step_size),
                        bn.DRAW_VARIANCE => string(hmcxml.draw_variance),
                        bn.AUTO_OPTIMIZE => string(hmcxml.auto_optimize),
                        bn.GRAD_CHECK_COUNT => string(hmcxml.check_grad),
                        bn.GRAD_TOLERANCE => string(hmcxml.grad_tolerance)])
    set_attributes(el, preconditioningAttrs(hmcxml.preconditioning))

    jg_el = new_child(el, bn.JOINT_GRADIENT)

    for i = 1:length(hmcxml.grads)
        grad = hmcxml.grads[i]
        make_xml(grad)
        if hmcxml.already_made[i]
            add_ref_el(jg_el, grad.el)
        else
            add_child(jg_el, grad.el)
        end
    end

    add_ref_el(el, get_hmc_parameter(hmcxml.param_provider))
    if hmcxml.transform != ""
        t_el = new_child(el, bn.TRANSFORM)
        set_attribute(t_el, bn.TYPE, hmcxml.transform)
    end

    if length(hmcxml.mask) > 0
        m_el = new_child(el, bn.MASK)
        add_child(m_el, make_parameter(value=hmcxml.mask))
    end


    hmcxml.el = el
    return el
end

function get_parameter(op::HMCOperatorXMLElement)
    return get_hmc_parameter(op.param_provider)
end

################################################################################
## Factor model loadings gradients
################################################################################


## loadings prior
mutable struct LoadingsGradientXMLElement <: MyXMLElement
    el::XMLOrNothing
    ifxml::IntegratedFactorsXMLElement

    function LoadingsGradientXMLElement(ifxml::IntegratedFactorsXMLElement)

        return new(nothing, ifxml)
    end

end

function make_xml(gxml::LoadingsGradientXMLElement)
    ifxml = gxml.ifxml
    el = make_loadings_gradient(ifxml.loadings, ifxml.loadings_prior)

    gxml.el = el
    return el
end

## integrated likelihood
mutable struct FactorLoadingsGradientXMLElement <: MyXMLElement
    el::XMLOrNothing
    ifxml::IntegratedFactorsXMLElement
    tdlxml::TraitLikelihoodXMLElement
    name::String


    function FactorLoadingsGradientXMLElement(
            ifxml::IntegratedFactorsXMLElement,
            tdlxml::TraitLikelihoodXMLElement)

        return new(nothing, ifxml, tdlxml, bn.INTEGRATED_FACTOR_GRADIENT)
    end

end

function make_xml(flgxml::FactorLoadingsGradientXMLElement)
    make_xml(flgxml.ifxml)
    make_xml(flgxml.tdlxml)

    el = new_element(flgxml.name)
    add_ref_el(el, flgxml.ifxml.el)
    add_ref_el(el, flgxml.tdlxml.el)

    flgxml.el = el
    return el
end



mutable struct ScaledMatrixGradient <: MyXMLElement
    el::XMLOrNothing
    original_gradient::MyXMLElement
    component::String

    function ScaledMatrixGradient(
            original_gradient::MyXMLElement,
            component::String)
        return new(nothing, original_gradient, component)
    end
end

function make_xml(smg::ScaledMatrixGradient)
    el = new_element(bn.SCALED_MATRIX_GRADIENT)
    set_attribute(el, bn.COMPONENT, smg.component)
    sub_el = make_xml(smg.original_gradient)
    add_child(el, sub_el)
    smg.el = el
    return el
end

function NormalizedLoadingsGradientXMLElement(
    ifxml::IntegratedFactorsXMLElement,
    tdlxml::TraitLikelihoodXMLElement)
    grad = FactorLoadingsGradientXMLElement(ifxml, tdlxml)
    return ScaledMatrixGradient(grad, bn.MATRIX)
end

function ScaleLoadingsGradientXMLElement(
    ifxml::IntegratedFactorsXMLElement,
    tdlxml::TraitLikelihoodXMLElement)
    grad = FactorLoadingsGradientXMLElement(ifxml, tdlxml)
    return ScaledMatrixGradient(grad, bn.SCALE)
end

## sampled likelihood
mutable struct SampledLoadingsGradient <: MyXMLElement
    el::XMLOrNothing
    factor_model::LatentFactorModelXMLElement

    function SampledLoadingsGradient(
            factor_model::LatentFactorModelXMLElement)
        return new(nothing, factor_model)
    end
end


function make_xml(slg::SampledLoadingsGradient)
    el = new_element(bn.SAMPLED_LOADINGS_GRADIENT)
    fac_el = make_xml(slg.factor_model)
    add_ref_el(el, fac_el)
    slg.el = el
    return el
end


################################################################################
## preconditioning
################################################################################


mutable struct PreconditioningSchedule
    delay::Int
    update_frequency::Int
    max_update::Int
    memory::Int
end

function preconditioningAttrs(ps::PreconditioningSchedule)
    return Dict{String, String}(
        bn.PRECONDITIONING_UPDATE_FREQUENCY => string(ps.update_frequency),
        bn.PRECONDITIONING_MAX_UPDATE => string(ps.max_update),
        bn.PRECONDITIONING_DELAY => string(ps.delay),
        bn.PRECONDITIONING_MEMORY => string(ps.memory)
    )
end

function PreconditioningSchedule()
    return PreconditioningSchedule(100, 100, 0, 100)
end

struct NoPreconditioning <: Preconditioner
end

function preconditioningAttrs(::NoPreconditioning)
    return Dict{String, String}()
end

mutable struct AdaptiveDiagonalPreconditioning <: Preconditioner
    schedule::PreconditioningSchedule
end

function AdaptiveDiagonalPreconditioning()
    return AdaptiveDiagonalPreconditioning(PreconditioningSchedule())
end

function preconditioningAttrs(adp::AdaptiveDiagonalPreconditioning)
    attrs = preconditioningAttrs(adp.schedule)
    attrs[bn.PRECONDITIONING] = bn.ADAPTIVE_DIAGONAL
    return attrs
end



