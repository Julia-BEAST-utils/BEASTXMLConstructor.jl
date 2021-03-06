mutable struct MultiplicativeScalePrior <: MyXMLElement
    mults::Parameter
    precs::MultiplicativeParameter
    mult_prior::MultivariateGammaLikelihood
    scale_prior::MyXMLElement
end

function make_xml(mps::MultiplicativeScalePrior)
    make_xml(mps.mults)
    make_xml(mps.precs)
    make_xml(mps.mult_prior)
    make_xml(mps.scale_prior)
    return nothing
end

function get_loggables(mps::MultiplicativeScalePrior)
    make_xml(mps)
    return [mps.mults.el]
end

function get_loadings_prior(msp::MultiplicativeScalePrior)
    return msp.scale_prior
end

function NormalGammaPrecisionOperatorXMLElement(
    multiplicativeLikelihood::MultiplicativeScalePrior)


    ggp = MultipilcativeGammaGibbsProvider(multiplicativeLikelihood.mults,
                multiplicativeLikelihood.scale_prior)

    return NormalGammaPrecisionOperatorXMLElement(nothing,
                ggp,
                multiplicativeLikelihood.mult_prior,
                Int[],
                1.0)
end

function set_shrinkage_mults!(msp::MultiplicativeScalePrior;
                        shapes::Vector{Float64} = Float64[],
                        scales::Vector{Float64} = Float64[],
                        set_mults::Bool = true)
    μ = shapes .* scales
    set_shrinkage_mults!(msp.mult_prior, shapes=shapes, scales=scales)
    if set_mults
        set_value(msp.mults, μ)
    end
    scale_mode = [prod(μ[1:i]) for i = 1:length(μ)]
    scale_mode = scale_mode.^(-0.5)
    return scale_mode
end

function get_priors(msp::MultiplicativeScalePrior)
    make_xml(msp)
    return [msp.mult_prior.el, msp.scale_prior.el]
end

function get_normal_prior(msp::MultiplicativeScalePrior)
    return get_normal_prior(msp.scale_prior)
end