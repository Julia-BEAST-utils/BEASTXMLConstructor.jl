mutable struct IndependentNormalDistributionModel <: MyXMLElement
    el::XMLOrNothing
    data::MyXMLElement
    mean::MyXMLElement
    precision::MyXMLElement
    make_data::Bool
    make_mean::Bool
    make_precision::Bool
    id::String

    function IndependentNormalDistributionModel(data::MyXMLElement,
                                                mean::MyXMLElement,
                                                precision::MyXMLElement,
                                                make_data::Bool,
                                                make_mean::Bool,
                                                make_precision::Bool,
                                                id::String)
        return new(nothing, data, mean, precision,
                   make_data, make_mean, make_precision, id)
    end

end

function check_parameter(x::MyXMLElement)
    return x, false
end

function check_parameter(x::Vector{Float64})
    return Parameter(x), true
end

function IndependentNormalDistributionModel(
            data::Union{MyXMLElement, Vector{Float64}}, id::String;
            mean::Union{MyXMLElement, Vector{Float64}} = zeros(length(data)),
            precision::Union{MyXMLElement, Vector{Float64}} = zeros(length(data)))

    data, make_data = check_parameter(data)
    mean, make_mean = check_parameter(mean)
    precision, make_precision = check_parameter(precision)

    return IndependentNormalDistributionModel(data, mean, precision,
                                              make_data,
                                              make_mean,
                                              make_precision, id)
end

function name(::IndependentNormalDistributionModel)
    return bn.INDEPENDENT_NORMAL_DISTRIBUTION
end

function get_loadings_prior(indnorm::IndependentNormalDistributionModel)
    return indnorm
end

function make_xml(indnorm::IndependentNormalDistributionModel)
    el = new_element(bn.INDEPENDENT_NORMAL_DISTRIBUTION)
    set_attribute(el, bn.ID, indnorm.id)
    data_el = new_child(el, bn.DATA)
    if indnorm.make_data
        add_child(data_el, make_xml(indnorm.data))
    else
        add_ref_el(data_el, indnorm.data)
    end

    mean_el = new_child(el, bn.MEAN)
    if indnorm.make_mean
        add_child(mean_el, make_xml(indnorm.mean))
    else
        add_ref_el(mean_el, indnorm.mean)
    end

    prec_el = new_child(el, bn.PRECISION)
    if indnorm.make_precision
        add_child(prec_el, make_xml(indnorm.precision))
    else
        add_ref_el(prec_el, indnorm.precision)
    end

    indnorm.el = el
    return el
end



