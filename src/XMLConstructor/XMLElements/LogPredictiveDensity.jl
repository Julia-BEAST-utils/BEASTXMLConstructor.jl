mutable struct FactorLogPredictiveDensity <: MyXMLElement
    fac_el::XMLOrNothing
    like_el::XMLOrNothing
    cl_el::XMLOrNothing # compound likelihood
    intfac::IntegratedFactorsXMLElement
    like::TraitLikelihoodXMLElement
    trait_ind::Int

    function FactorLogPredictiveDensity(
                intfac::IntegratedFactorsXMLElement,
                like::TraitLikelihoodXMLElement;
                trait_ind::Int = 2)

        return new(nothing, nothing, nothing, intfac, like, trait_ind)
    end
end

function make_xml(lpd::FactorLogPredictiveDensity)

    intfac = copy(lpd.intfac)
    intfac.tree_trait_ind = lpd.trait_ind # TODO: assumes only two traits
    like = copy(lpd.like)
    trait_name = intfac.treeModel.node_traits[intfac.tree_trait_ind]
    intfac.id = "$(trait_name).factorLikelihood"
    intfac.standardize_traits = false
    fac_el = make_xml(intfac, reference_precision = true)
    lpd.fac_el = fac_el


    like.extension_el = intfac
    lpd.like_el = make_xml(like)
    set_id!(lpd.like_el, "$(trait_name).treeLikelihood")

    cl_el = new_element(bn.LIKELIHOOD)
    set_id!(cl_el, "$trait_name.likelihood")
    add_ref_el(cl_el, lpd.fac_el)
    add_ref_el(cl_el, lpd.like_el)
    lpd.cl_el = cl_el

end

function get_loggables(lpd::FactorLogPredictiveDensity)
    make_xml(lpd)
    return lpd.cl_el
end
