abstract type GammaGibbsProvider <: MyXMLElement end

mutable struct NormalGammaPrecisionOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    ggp::GammaGibbsProvider
    prior_provider::MyXMLElement
    indices::Vector{Int}
    weight::Float64
end

function NormalGammaPrecisionOperatorXMLElement(ggp::GammaGibbsProvider,
                    prior::T) where T <: MyXMLElement
    return NormalGammaPrecisionOperatorXMLElement(nothing, ggp, prior, Int[], 1.0)
end


function NormalGammaPrecisionOperatorXMLElement(
            ext_xml::ModelExtensionXMLElement,
            td_xml::TraitLikelihoodXMLElement)
    return NormalGammaPrecisionOperatorXMLElement(nothing,
                NormalExtentsionGibbsProvider(ext_xml, td_xml),
                ext_xml, Int[], 1.0)
end

# function NormalGammaPrecisionOperatorXMLElement(msl::MatrixShrinkageLikelihoods,
#                                                 row::Int)
#     gpp = MultipilcativeGammaGibbsProvider(msl, row)
#     return NormalGammaPrecisionOperatorXMLElement(nothing, gpp, gpp, 1.0)
# end

function get_parameter(ngpo::NormalGammaPrecisionOperatorXMLElement)
    return get_parameter(ngpo.ggp)
end

mutable struct NormalExtentsionGibbsProvider <: GammaGibbsProvider
    el::XMLOrNothing
    ext::ModelExtensionXMLElement
    tdl::TraitLikelihoodXMLElement

    function NormalExtentsionGibbsProvider(ext::ModelExtensionXMLElement,
                                           tdl::TraitLikelihoodXMLElement)
        return new(nothing, ext, tdl)
    end
end

function get_parameter(negp::NormalExtentsionGibbsProvider)
    return get_precision(negp.ext)
end

function make_xml(negp::NormalExtentsionGibbsProvider)
    el = new_element(bn.NORMAL_EXTENSION)

    set_attribute(el, bn.TREE_TRAIT_NAME, negp.tdl.attrs[bn.TRAIT_NAME])
    add_ref_el(el, negp.ext)
    add_ref_el(el, negp.tdl)

    negp.el = el
    return el
end

mutable struct MultipilcativeGammaGibbsProvider <: GammaGibbsProvider
    el::XMLOrNothing
    param::MyXMLElement
    likelihood::MyXMLElement

    function MultipilcativeGammaGibbsProvider(param::MyXMLElement,
                                              likelihood::MyXMLElement)
        return new(nothing, param, likelihood)
    end
end

# function multiplicative_gamma_gibbs_proivders(msl::MatrixShrinkageLikelihoods)
#     k = get_fac_dim(msl)
#     mggs = Vector{MultipilcativeGammaGibbsProvider}(undef, k)
#     for i = 1:k
#         mggs[i] = MultipilcativeGammaGibbsProvider(msl, i)
#     end
#     return mggs
# end

function get_parameter(mggp::MultipilcativeGammaGibbsProvider)
    return mggp.param
end

function make_xml(mgg::MultipilcativeGammaGibbsProvider)
    make_xml(mgg.param)
    make_xml(mgg.likelihood)

    el = new_element(bn.MULTIPLICATIVE_GAMMA_GIBBS)
    add_ref_el(el, mgg.param.el)
    add_ref_el(el, mgg.likelihood.el)

    mgg.el = el
end

# function get_precision_prior(mgg::MultipilcativeGammaGibbsProvider)
#     make_xml(mgg.msl)
#     return mgg.msl.global_prior_els[mgg.row]
# end




function make_xml(ngpxml::NormalGammaPrecisionOperatorXMLElement)
    make_xml(ngpxml.ggp)
    make_xml(ngpxml.prior_provider)

    el = new_element(bn.NGP_OPERATOR)
    set_attribute(el, bn.WEIGHT, ngpxml.weight)
    if length(ngpxml.indices) > 0
        set_attribute(el, bn.INDICES, join(ngpxml.indices, ' '))
    end
    prior_el = new_child(el, bn.PRIOR)
    add_ref_el(prior_el, get_precision_prior(ngpxml.prior_provider))

    add_child(el, ngpxml.ggp.el)

    ngpxml.el = el
    return el
end

function set_indices!(ngp::NormalGammaPrecisionOperatorXMLElement,
                     inds::Vector{Int})
    ngp.indices = inds
end
