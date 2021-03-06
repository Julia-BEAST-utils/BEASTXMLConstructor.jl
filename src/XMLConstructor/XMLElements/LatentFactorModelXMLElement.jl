mutable struct LatentFactorModelXMLElement <: ModelExtensionXMLElement
    el::XMLOrNothing
    precision_prior_el::XMLOrNothing
    loadings::MyXMLElement
    loadings_prior::MyXMLElement
    precision::MyXMLElement
    precision_scale::Float64
    precision_shape::Float64
    treeModel::TreeModelXMLElement
    tdl::TraitLikelihoodXMLElement
    tree_trait_ind::Int
    trait_name::String
    standardize_traits::Bool
    parameters_already_made::Bool
    id::String

    function LatentFactorModelXMLElement(treeModel_xml::TreeModelXMLElement,
                tl_xml::TraitLikelihoodXMLElement, k::Int)


        p = treeModel_xml.trait_dims[1]

        L = default_loadings(k, p)
        L_param = MatrixParameter(L, bn.DEFAULT_LOADINGS_ID)

        trait_name = treeModel_xml.node_traits[1]
        precision = Parameter(ones(p), bn.FACTOR_PRECISION)
        precision_scale = 1.0
        precision_shape = 1.0
        tree_trait_ind = 1
        return new(nothing, nothing,
                L_param,
                NormalDistributionLikelihood(
                                        L_param,
                                        bn.DEFAULT_LOADINGS_ID * ".prior"),
                precision, precision_scale,
                precision_shape,
                treeModel_xml,
                tl_xml,
                tree_trait_ind, trait_name, true, false,
                bn.LATENT_FACTOR_MODEL)
    end
end

# function default_loadings(k::Int, p::Int)
#     @assert k <= p
#     L = zeros(k, p)
#     for i = 1:k
#         for j = i:p
#             L[i, j] = 1.0
#         end
#     end
#     # L = randn(k, p)
#     return L
# end
function name(::LatentFactorModelXMLElement)
    return bn.LATENT_FACTOR_MODEL
end

function make_xml(lfxml::LatentFactorModelXMLElement)


    # set_attribute(lfxml.loadings_el, bn.ID, bn.LOADINGS)

    el = new_element(bn.LATENT_FACTOR_MODEL)
    set_attribute(el, bn.ID, bn.LATENT_FACTOR_MODEL)
    set_attributes(el, [bn.FACTOR_NUMBER => string(size(lfxml.loadings, 1)),
                        bn.TRAIT_NAME => lfxml.trait_name,
                        bn.SCALE_DATA => string(lfxml.standardize_traits)])
    f_el = new_child(el, bn.FACTORS)

    make_xml(lfxml.tdl)
    tp_el = find_element(find_element(lfxml.tdl.el, bn.TRAIT_PARAMETER), bn.PARAMETER)
    add_ref_el(f_el, tp_el)

    data_el = new_child(el, bn.DATA)
    dtt_el = new_child(data_el, bn.DATA_FROM_TIPS)
    set_attribute(dtt_el, bn.TRAIT_NAME, lfxml.trait_name)

    make_xml(lfxml.treeModel)
    add_ref_el(dtt_el, lfxml.treeModel.el)
    tp_el = new_child(dtt_el, bn.TRAIT_PARAMETER)
    add_ref_el(tp_el, bn.PARAMETER, lfxml.treeModel.param_names[1])

    add_ref_el(el, lfxml.treeModel.el)
    l_el = new_child(el,bn.LOADINGS)
    add_ref_el(l_el, lfxml.loadings)

    rp_el = new_child(el, bn.ROW_PRECISION)
    dm_el = new_child(rp_el, bn.DIAGONAL_MATRIX)
    add_parameter(dm_el, value=[1], dim=size(lfxml.loadings,2))

    cp_el = new_child(el, bn.COL_PRECISION)
    dm_el = new_child(cp_el, bn.DIAGONAL_MATRIX)
    if lfxml.parameters_already_made
        add_ref_el(dm_el, bn.PARAMETER, bn.FACTOR_PRECISION)
    else
        prec_el = make_xml(lfxml.precision)
        add_child(dm_el, prec_el)
    end



    lfxml.el = el

    lfxml.precision_prior_el = make_precision_prior(lfxml)

    # lfxml.loadings_prior_el = make_loadings_prior(lfxml)
    return el
end

function make_precision_prior(lfxml::LatentFactorModelXMLElement)
    el = new_element(bn.GAMMA_PRIOR)
    set_attributes(el, [bn.ID => "$(bn.FACTOR_PRECISION).prior",
                        bn.SCALE => string(lfxml.precision_scale),
                        bn.SHAPE => string(lfxml.precision_shape)])
    add_ref_el(el, lfxml.el[bn.COL_PRECISION][1][bn.DIAGONAL_MATRIX][1][bn.PARAMETER][1])
    return el
end

function get_loadings(lfm::LatentFactorModelXMLElement)
    return lfm.loadings
end

function get_precision(lfm::LatentFactorModelXMLElement)
    return lfm.precision
end


# function make_lf_loadings(L::AbstractArray{Float64, 2})
#     el = new_element(bn.MATRIX_PARAMETER)
#     set_attribute(el, bn.ID, bn.DEFAULT_LOADINGS_ID)
#     n, p = size(L)
#     for i = 1:n
#         add_parameter(el, id = "$(bn.DEFAULT_LOADINGS_ID)$i",
#                         value = L[i, :])
#     end
#     return el
# end

# mutable struct NormalDistributionLikelihood <: MyXMLElement
#     el::XMLOrNothing
#     parameter::MyXMLElement
#     μ::Float64
#     σ::Float64

#     function NormalDistributionLikelihood(parameter::MyXMLElement)
#         return new(nothing, parameter, 0.0, 1.0)
#     end
# end

# function make_xml(ndl::NormalDistributionLikelihood)
#     el = new_element(bn.DISTRIBUTION_LIKELIHOOD)
#     set_attribute(el, bn.ID, "$(get_id(ndl.parameter)).prior")

#     d_el = new_child(el, bn.DATA)
#     add_ref_el(d_el, ndl.parameter)
#     dist_el = new_child(el, bn.DISTRIBUTION)

#     ndm_el = new_child(dist_el, bn.NORMAL_DISTRIBUTION_MODEL)

#     m_el = new_child(ndm_el, bn.MEAN)
#     add_parameter(m_el, value = [μ])
#     std_el = new_child(ndm_el, bn.STDEV)
#     add_parameter(std_el, value =[σ], lower=0.0)

#     return el
# end


# function make_loadings_prior(lfxml::LatentFactorModelXMLElement)
#     el = new_element(bn.DISTRIBUTION_LIKELIHOOD)
#     set_attribute(el, bn.ID, "$(bn.DEFAULT_LOADINGS_ID).prior")

#     d_el = new_child(el, bn.DATA)
#     add_ref_el(d_el, lfxml.loadings_el)
#     dist_el = new_child(el, bn.DISTRIBUTION)

#     ndm_el = new_child(dist_el, bn.NORMAL_DISTRIBUTION_MODEL)

#     m_el = new_child(ndm_el, bn.MEAN)
#     add_parameter(m_el, value = [0.0])
#     std_el = new_child(ndm_el, bn.STDEV)
#     add_parameter(std_el, value =[1.0], lower=0.0)

#     return el
# end

function get_priors(xml::LatentFactorModelXMLElement)
    make_xml(xml)
    return [get_priors(xml.loadings_prior); xml.precision_prior_el]
end

function get_loggables(xml::LatentFactorModelXMLElement)
    make_xml(xml)
    return [get_loggables(xml.loadings);
            [xml.el[bn.COL_PRECISION][1][bn.DIAGONAL_MATRIX][1][bn.PARAMETER][1]]]
end

#
# function get_precision_prior(xml::IntegratedFactorsXMLElement)
#     return xml.precision_prior_el
# end
