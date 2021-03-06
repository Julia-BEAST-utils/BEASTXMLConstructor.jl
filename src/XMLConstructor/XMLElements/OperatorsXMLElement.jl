mutable struct OperatorsXMLElement <: MyXMLElement
    el::XMLOrNothing
    els::Vector{<:MyXMLElement}

    OperatorsXMLElement(el::MyXMLElement) = new(nothing, [el])
    OperatorsXMLElement(els::Vector{<:MyXMLElement}) = new(nothing, els)

end

function make_xml(os_el::OperatorsXMLElement)
    for el in os_el.els
        make_xml(el)
    end

    el = new_element(bn.OPERATORS)
    set_attribute(el, bn.ID, bn.OPERATORS)
    for o_el in os_el.els
        if typeof(o_el) <: CompoundOperatorXMLElement
            for xml_el in o_el.els
                add_child(el, xml_el)
            end
        else
            add_child(el, o_el.el)
        end
    end
    os_el.el = el
    return el
end

function merge_operators!(ops1::OperatorsXMLElement, ops2::OperatorsXMLElement)
    ops1.els = [ops1.els; ops2.els]
    return ops1
end

mutable struct FireParameterOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    param_provider::MyXMLElement
    value::Array{T} where T <: Real
    weight::Float64

    FireParameterOperatorXMLElement(my_xml::MyXMLElement) =
            new(nothing, my_xml, zeros(0), 1.0)
end

function make_xml(fpo_el::FireParameterOperatorXMLElement)

    make_xml(fpo_el.param_provider)

    el = new_element(bn.FIRE_PARAMETER_CHANGED)
    set_attribute(el, bn.WEIGHT, fpo_el.weight)

    if length(fpo_el.value) > 0
        value_string = join(fpo_el.value, ' ')
        set_attribute(el, bn.VALUE, value_string)
    end

    param_el = get_matrix_parameter(fpo_el.param_provider)
    add_ref_el(el, param_el)

    fpo_el.el = el

    return el
end

################################################################################
## scale operator
################################################################################

mutable struct ScaleOperator <: OperatorXMLElement
    el::XMLOrNothing
    param::MyXMLElement
    scale_factor::Float64
    inds::Vector{Int}
    weight::Float64

    function ScaleOperator(param::MyXMLElement, scale_factor::Float64,
                           weight::Float64)
        return new(nothing, param, scale_factor, Int[], weight)
    end
end

function ScaleOperator(param::MyXMLElement)
    return ScaleOperator(param, 0.5, 1.0)
end

function make_xml(so::ScaleOperator)
    el = new_element(bn.SCALE_OPERATOR)
    set_attribute(el, bn.SCALE_FACTOR, so.scale_factor)
    set_attribute(el, bn.WEIGHT, so.weight)
    if length(so.inds) == 0
        add_ref_el(el, so.param)
    else
        mask = zeros(length(so.param))
        mask[so.inds] .= 1.0
        mask_param = MaskedParameter(so.param, mask)
        add_child(el, make_xml(mask_param))
    end

    so.el = el
    return el
end


################################################################################
## shrinkage operators
################################################################################





abstract type CompoundOperatorXMLElement <: OperatorXMLElement end

mutable struct ShrinkageScaleOperators <: CompoundOperatorXMLElement
    els::Vector{XMLOrNothing}
    msl::MatrixShrinkageLikelihoods
    ifxml::IntegratedFactorsXMLElement
    global_weights::Vector{Float64}
    bb_weights::Vector{Float64}
    fix_first::Bool
    fix_globals::Bool

    function ShrinkageScaleOperators(msl::MatrixShrinkageLikelihoods,
                                    ifxml::IntegratedFactorsXMLElement)
        k = get_fac_dim(msl)
        return new(xml_vec(2 * k),
                    msl,
                    ifxml,
                    ones(k),
                    ones(k),
                    false,
                    false)
    end
end

function make_scale_operator(param::XMLElement,
                            scale_factor::Float64,
                            weight::Float64)
    return make_xml(ScaleOperator(param, scale_factor, weight))
end

function make_bb_operator(bb::XMLElement, weight::Float64)
    el = new_element(bn.BAYESIAN_BRIDGE_GIBBS_OP)
    set_attribute(el, bn.WEIGHT, weight)
    add_ref_el(el, bb)
    return el
end

function make_xml(xml::ShrinkageScaleOperators)
    msl = xml.msl
    make_xml(xml.ifxml)
    make_xml(msl)
    k = get_fac_dim(msl)
    if xml.fix_globals
        xml.els[1:k] .= nothing
    else
        if xml.fix_first
            xml.els[1] = nothing
        else
            op = NormalGammaPrecisionOperatorXMLElement(xml.msl, 1)
            make_xml(op)
            xml.els[1] = op.el
            # xml.els[1] = make_scale_operator(msl.gp_els[1], xml.scale_factors[1],
            #                                 xml.scale_weights[1])
        end

        for i = 2:k
            gibbs_provider = MultipilcativeGammaGibbsProvider(xml.msl, i)
            prior_provider = gibbs_provider
            op = NormalGammaPrecisionOperatorXMLElement(gibbs_provider,
                                                        prior_provider)
            make_xml(op)
            xml.els[i] = op.el
            # xml.els[i] = make_scale_operator(msl.mult_els[i - 1],
            #                                 xml.scale_factors[i],
            #                                 xml.scale_weights[i])
        end
    end

    for i = 1:k
        xml.els[i + k] = make_bb_operator(msl.bb_els[i], xml.bb_weights[i])
    end

end


################################################################################
## joint operator
################################################################################

mutable struct JointOperator <: OperatorXMLElement
    el::XMLOrNothing
    operators::Vector{<:OperatorXMLElement}
    weight::Float64

    function JointOperator(operators::Vector{<:OperatorXMLElement},
                           weight::Float64)
        return new(nothing, operators, weight)
    end
end

function make_xml(jo::JointOperator)
    el = new_element(bn.JOINT_OPERATOR)
    for op in jo.operators
        op_el = make_xml(op)
        add_child(el, op_el)
    end
    set_attribute(el, bn.WEIGHT, jo.weight)
    jo.el = el
    return el
end

################################################################################
## rejection operator
################################################################################

mutable struct RejectionOperator <: OperatorXMLElement
    el::XMLOrNothing
    operator::OperatorXMLElement
    condition::String

    function RejectionOperator(op::OperatorXMLElement, condition::String)
        return new(nothing, op, condition)
    end
end

function name(::RejectionOperator)
    return bn.REJECTION_OPERATOR
end

function make_xml(ro::RejectionOperator)
    el = new_element(name(ro))
    set_attribute(el, bn.CONDITION, ro.condition)
    add_child(el, make_xml(ro.operator))
    ro.el = el
    return el
end