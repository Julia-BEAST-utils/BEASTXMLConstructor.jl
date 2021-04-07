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
                trait_name::String = bn.DEFAULT_TRAIT_NAME)

        name_ind = findfirst(isequal(trait_name), intfac.treeModel.node_traits)

        intfac = copy(intfac)
        intfac.tree_trait_ind = name_ind
        like = copy(like)
        like.attrs[bn.ID] = "$(trait_name).treeLikelihood"
        intfac.id = "$(trait_name).factorLikelihood"
        intfac.standardize_traits = false


        like.extension_el = intfac

        return new(nothing, nothing, nothing, intfac, like, name_ind)
    end
end

function make_xml(lpd::FactorLogPredictiveDensity)

    @unpack intfac, like = lpd
    fac_el = make_xml(intfac, reference_precision = true)
    lpd.fac_el = fac_el
    trait_name = get_trait_name(intfac)


    like.extension_el = intfac
    lpd.like_el = make_xml(like)

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
