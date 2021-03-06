mutable struct MultivariateGammaLikelihood <: MyXMLElement
    el::XMLOrNothing
    data::Parameter
    shapes::Vector{Float64}
    scales::Vector{Float64}
    id::String

    function MultivariateGammaLikelihood(data::Parameter,
                                         shapes::Vector{Float64},
                                         scales::Vector{Float64},
                                         id::String)
        @assert length(data.val) == length(shapes) == length(scales)
        return new(nothing, data, shapes, scales, id)
    end
end

function make_xml(mgl::MultivariateGammaLikelihood)
    make_xml(mgl.data)

    el = new_element(bn.MULTIVARIATE_GAMMA_LIKELIHOOD)
    set_id!(el, mgl.id)

    data_el = new_child(el, bn.DATA)
    add_ref_el(data_el, mgl.data.el)

    shape_el = new_child(el, bn.SHAPE)
    add_parameter(shape_el, value = mgl.shapes, lower = 0.0)

    scale_el = new_child(el, bn.SCALE)
    add_parameter(scale_el, value = mgl.scales, lower = 0.0)

    mgl.el = el
    return el
end

function get_precision_prior(mgl::MultivariateGammaLikelihood)
    make_xml(mgl)
    return mgl.el
end

function set_shrinkage_mults!(mgl::MultivariateGammaLikelihood;
                            shapes::Vector{Float64} = Float64[],
                            scales::Vector{Float64} = Float64[])

    n = length(mgl.shapes)
    set_values(shapes, mgl.shapes)
    set_values(scales, mgl.scales)
end


function set_values(src::Vector{Float64}, dest::Vector{Float64})
    @show src
    @show dest
    if length(src) > 0
        if length(src) == length(dest) - 1
            dest[2:end] .= src
        else
            dest .= src
        end
    end
end