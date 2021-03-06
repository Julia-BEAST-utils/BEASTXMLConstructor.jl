const DEFAULT_SHRINKAGE_SHAPE = 100.0
const DEFAULT_SHRINKAGE_SCALE = 1.0
const DEFAULT_STARTING_SHAPE = 2.0


mutable struct MatrixShrinkageLikelihoods <: MyXMLElement
    mult_els::Vector{XMLOrNothing}
    gp_els::Vector{XMLOrNothing}
    local_scale_els::Vector{XMLOrNothing}
    global_scale_els::Vector{XMLOrNothing}
    bb_els::Vector{XMLOrNothing}
    ms_el::XMLOrNothing
    global_prior_els::Vector{XMLOrNothing}
    loadings::MatrixParameter
    shapes::Vector{Float64}
    scales::Vector{Float64}
    local_scales::Matrix{Float64}

    function MatrixShrinkageLikelihoods(loadings::MatrixParameter)

        k, p = size(loadings.mat)
        shapes = fill(DEFAULT_SHRINKAGE_SHAPE, k)
        shapes[1] = DEFAULT_STARTING_SHAPE
        scales = fill(DEFAULT_SHRINKAGE_SCALE, k)

        return new(xml_vec(k - 1),
                    xml_vec(k),
                    xml_vec(k),
                    xml_vec(k),
                    xml_vec(k),
                    nothing,
                    xml_vec(k),
                    loadings,
                    shapes,
                    scales,
                    ones(k, p)
                    )
    end
end

function get_fac_dim(msl::MatrixShrinkageLikelihoods)
    return length(msl.gp_els)
end


function make_xml(msl::MatrixShrinkageLikelihoods)

    loadings = make_xml(msl.loadings)
    load_rows = loadings[bn.PARAMETER]
    k = length(load_rows)
    p = length(split(attribute(load_rows[1], bn.VALUE)))

    l_ids = Vector{String}(undef, k)
    for i = 1:k
        l_ids[i] = attribute(load_rows[i], bn.ID)
    end

    for i = 2:k
        msl.mult_els[i - 1] =
                make_parameter(id="rowMult$i",
                                value=[msl.shapes[i] * msl.scales[i]],
                                lower=0.0)
    end

    for i = 1:k
        msl.local_scale_els[i] = make_parameter(id="localScale$i",
                                    value=vec(msl.local_scales[i, :]),
                                    lower=0.0)
    end

    msl.gp_els[1] = make_parameter(id="globalPrecision1",
                                    value=[msl.shapes[1] * msl.scales[1]],
                                    lower=0.0)

    for i = 2:k
        prod_param = new_element(bn.PRODUCT_PARAMETER)
        set_id!(prod_param, "globalPrecision$i")
        add_ref_el(prod_param, msl.gp_els[i - 1])
        add_ref_el(prod_param, msl.mult_els[i - 1])
        msl.gp_els[i] = prod_param
    end

    for i = 1:k
        trans_param = new_element(bn.TRANSFORMED_PARAMETER)
        set_id!(trans_param, "globalScale$i")
        set_attribute(trans_param, bn.TYPE, bn.POWER)
        add_ref_el(trans_param, msl.gp_els[i])
        pt_el = new_child(trans_param, bn.POWER_TRANSFORM)
        set_attribute(pt_el, bn.POWER, "-0.5")
        t_el = new_child(pt_el, bn.TRANSFORM)
        set_attribute(t_el, bn.TYPE, bn.POWER)
        msl.global_scale_els[i] = trans_param
    end

    for i = 1:k
        bb_el = new_element(bn.BAYESIAN_BRIDGE)
        set_id!(bb_el, "bb$i")
        add_ref_el(bb_el, load_rows[i])
        gs_el = new_child(bb_el, bn.GLOBAL_SCALE)
        add_ref_el(gs_el, msl.global_scale_els[i])
        ls_el = new_child(bb_el, bn.LOCAL_SCALE)
        add_ref_el(ls_el, msl.local_scale_els[i])
        exp_el = new_child(bb_el, bn.EXPONENT)
        add_parameter(exp_el, value=[0.5])
        msl.bb_els[i] = bb_el
    end

    ms_el = new_element(bn.MATRIX_SHRINKAGE_LIKELIHOOD)
    set_id!(ms_el, "$(get_id(loadings)).prior")
    add_ref_el(ms_el, loadings)
    rows_el = new_child(ms_el, bn.ROW_PRIORS)
    for i = 1:k
        add_ref_el(rows_el, msl.bb_els[i])
    end
    msl.ms_el = ms_el

    for i = 1:k

        gp_prior = new_element(bn.GAMMA_PRIOR)
        if i == 1
            set_id!(gp_prior, "globalPrecisionPrior1")
        else
            set_id!(gp_prior, "rowMultPrior$i")
        end
        set_attribute(gp_prior, bn.SHAPE, msl.shapes[i])
        set_attribute(gp_prior, bn.SCALE, msl.scales[i])
        if i == 1
            add_ref_el(gp_prior, msl.gp_els[1])
        else
            add_ref_el(gp_prior, msl.mult_els[i - 1])
        end

        msl.global_prior_els[i] = gp_prior
    end




    all_els = [msl.mult_els; msl.gp_els; msl.local_scale_els;
                msl.global_scale_els; msl.bb_els; msl.ms_el; msl.global_prior_els]

    return all_els

end

function get_priors(xml::MatrixShrinkageLikelihoods)
    return XMLElement[xml.ms_el; xml.global_prior_els]
end

function get_loggables(xml::MatrixShrinkageLikelihoods)
    return [xml.global_scale_els; xml.mult_els; xml.local_scale_els]
end

function get_precision_prior(xml::MatrixShrinkageLikelihoods, ind::Int; make_xml::Bool = true)
    if make_xml
        make_xml(xml)
    end
    return xml.global_prior_els[ind]
end

function make_loadings_gradient(::MyXMLElement, msl::MatrixShrinkageLikelihoods)
    make_xml(msl)
    return msl.ms_el
end

function get_normal_prior(msl::MatrixShrinkageLikelihoods)
    make_xml(msl)
    return msl.ms_el
end