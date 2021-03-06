mutable struct WishartStatXMLElement <: MyXMLElement
    el::XMLOrNothing
    traitLikelihood_el::TraitLikelihoodXMLElement
    rm_el::Union{Nothing, RepeatedMeasuresXMLElement}
    mbd_el::Union{Nothing, MBDXMLElement}

    function WishartStatXMLElement(like_el::TraitLikelihoodXMLElement,
                mbd_el::MBDXMLElement)
        return new(nothing, like_el, nothing, mbd_el)
    end

    function WishartStatXMLElement(like_el::TraitLikelihoodXMLElement,
                rm_el::RepeatedMeasuresXMLElement)
        return new(nothing, like_el, rm_el, nothing)
    end
end

mutable struct PrecisionGibbsOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    wishart_stats::WishartStatXMLElement
    weight::Float64

    function PrecisionGibbsOperatorXMLElement(
            traitLikelihood_el::TraitLikelihoodXMLElement,
            mbd_el::MBDXMLElement)

        wishart_stats = WishartStatXMLElement(traitLikelihood_el, mbd_el)
        return new(nothing, wishart_stats, 1.0)
    end

    function PrecisionGibbsOperatorXMLElement(
        traitLikelihood_el::TraitLikelihoodXMLElement,
        rm_el::RepeatedMeasuresXMLElement)

        wishart_stats = WishartStatXMLElement(traitLikelihood_el, rm_el)
        return new(nothing, wishart_stats, 1.0)
    end

end

mutable struct CompoundPrecisionOperatorXMLElement <: OperatorXMLElement
    el::XMLOrNothing
    diff_el::PrecisionGibbsOperatorXMLElement
    res_el::PrecisionGibbsOperatorXMLElement
    weight::Float64

    function CompoundPrecisionOperatorXMLElement(
            diff_el::PrecisionGibbsOperatorXMLElement,
            res_el::PrecisionGibbsOperatorXMLElement)

        diff_ws = diff_el.wishart_stats
        res_ws = res_el.wishart_stats
        @assert !isnothing(diff_ws.mbd_el) && isnothing(diff_ws.rm_el)
        @assert isnothing(res_ws.mbd_el) && !isnothing(res_ws.rm_el)

        return new(nothing, diff_el, res_el, 1.0)
    end

    function CompoundPrecisionOperatorXMLElement(
                tl_el::TraitLikelihoodXMLElement,
                mbd_el::MBDXMLElement,
                rm_el::RepeatedMeasuresXMLElement)

        diff_el = PrecisionGibbsOperatorXMLElement(tl_el, mbd_el)
        res_el = PrecisionGibbsOperatorXMLElement(tl_el, rm_el)
        return CompoundPrecisionOperatorXMLElement(diff_el, res_el)
    end

end

function make_xml(cpo_el::CompoundPrecisionOperatorXMLElement)
    el = new_element(bn.COMPOUND_PRECISION_OPERATOR)
    set_attribute(el, bn.WEIGHT, cpo_el.weight)

    do_el = new_child(el, bn.DIFFUSION_OPERATOR)
    make_xml(cpo_el.diff_el)
    add_child(do_el, cpo_el.diff_el.el)

    ro_el = new_child(el, bn.RESIDUAL_OPERATOR)
    make_xml(cpo_el.res_el)
    add_child(ro_el, cpo_el.res_el.el)

    cpo_el.el = el
    return el
end




function make_xml(ws_el::WishartStatXMLElement)
    make_xml(ws_el.traitLikelihood_el)

    if isnothing(ws_el.rm_el) && !isnothing(ws_el.mbd_el)
        el = new_element(bn.WISHART_STATISTICS)
        return finish_Wishart_el(ws_el, el)
    elseif !isnothing(ws_el.rm_el) && isnothing(ws_el.mbd_el)
        make_xml(ws_el.rm_el)
        el = new_element(bn.REPEATED_MEASURES_WISHART_STATISTICS)
        add_ref_el(el, ws_el.rm_el.el)
        return finish_Wishart_el(ws_el, el)
    else
        error("Not implemented")
    end
end

function finish_Wishart_el(ws_el::WishartStatXMLElement, el::XMLElement)
    add_ref_el(el, ws_el.traitLikelihood_el.el)
    set_attribute(el, bn.TRAIT_NAME, ws_el.traitLikelihood_el.attrs[bn.TRAIT_NAME])
    ws_el.el = el
    return el
end

function make_xml(pgo_el::PrecisionGibbsOperatorXMLElement)
    el = new_element(bn.PRECISION_GIBBS_OPERATOR)
    set_attribute(el, bn.WEIGHT, pgo_el.weight)
    ws_el = pgo_el.wishart_stats
    make_xml(ws_el)
    add_child(el, ws_el.el)
    if isnothing(ws_el.rm_el) && !isnothing(ws_el.mbd_el)
        make_xml(ws_el.mbd_el)
        add_ref_el(el, ws_el.mbd_el.precision_prior.el)
    elseif !isnothing(ws_el.rm_el) && isnothing(ws_el.mbd_el)
        make_xml(ws_el.rm_el)
        add_ref_el(el, ws_el.rm_el.precision_prior.el)
    else
        error("Not implemented")
    end
    pgo_el.el = el
    return el
end
