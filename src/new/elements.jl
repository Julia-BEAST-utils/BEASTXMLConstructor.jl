
################################################################################
## FactorModel
################################################################################

function integratedFactorModelXML(;
        loadings::GeneralizedXMLElement,
        precision::GeneralizedXMLElement,
        treeModel::GeneralizedXMLElement,
        traitParameter::GeneralizedXMLElement,
        trait_name::String,
        id::Union{String, Nothing} = nothing,
        standardize::Bool = true
        )
    children = [PassthroughXMLElement("loadings", loadings),
                PassthroughXMLElement("precision", precision),
                treeModel,
                PassthroughXMLElement("traitParameter", traitParameter)]
    attributes = [bn.TRAIT_NAME => trait_name, bn.STANDARDIZE => "$standardize"]
    return GeneralizedXMLElement("integratedFactorModel", id, attributes, children, nothing)
end

function find_loadings(xml::GeneralizedXMLElement)
    @assert xml.name == "integratedFactorModel"
    return find_element(xml, name = "loadings", passthrough = true)
end

function find_factor_precision(xml::GeneralizedXMLElement)
    @assert xml.name == "integratedFactorModel"
    return find_element(xml, name="precision", passthrough = true)
end

################################################################################
## repeatedMeasures
################################################################################

"""
<repeatedMeasuresModel id="traits.residualModel" traitName="traits" forceFullPrecision="true">
    <treeModel idref="treeModel"/>
    <traitParameter>
        <parameter id="leaf.traits"/>
    </traitParameter>
    <samplingPrecision>

    </samplingPrecision>
</repeatedMeasuresModel>
"""
function repeatedMeasuresModelXML(;
            tree_model::GeneralizedXMLElement,
            precision::GeneralizedXMLElement,
            trait_name::String,
            id::String = "$trait_name.residualModel",
            standardize::Bool = false
            )
    trait_parameter = find_trait_parameter(tree_model, trait_name)

    return GeneralizedXMLElement("repeatedMeasuresModel", id = id,
            children = [
                    tree_model,
                    PassthroughXMLElement("traitParameter", trait_parameter),
                    PassthroughXMLElement("samplingPrecision", precision)
                    ],
            attributes = [bn.TRAIT_NAME => trait_name,
                          bn.STANDARDIZE => standardize,
                          bn.FORCE_FULL_PRECISION => true]
        )
end

################################################################################
## Data
################################################################################

function taxaXML(taxa::Vector{String},
        data::DataPairs;
        dates::Union{Nothing, Vector{Float64}} = nothing)

    n = length(taxa)
    p = length(data)

    children = Vector{GeneralizedXMLElement}(undef, n)
    for i = 1:n
        date = isnothing(dates) ? nothing : dates[i]
        traits = [data[j][1] => data[j][2][i, :] for j = 1:p]
        children[i] = taxonXML(taxa[i], traits, date)
    end
    return GeneralizedXMLElement("taxa", children = children)
end

function taxonXML(taxon::String,
        traits::Vector{<:Pair{String, <:AbstractVector{Float64}}},
        date::Union{Float64, Nothing})
    children = [traitAttr(trait) for trait in traits]
    if !isnothing(date)
        error("not yet implemented")
    end
    return GeneralizedXMLElement("taxon", id = taxon, children = children)
end

function traitAttr(attr::Pair{String, <:AbstractArray{Float64}})
    values = attr[2]
    n = length(values)
    formatted_values = fill("", n)
    for i = 1:n
        if isnan(values[i])
            formatted_values[i] = "NA"
        else
            formatted_values[i] = string(values[i])
        end
    end
    return GeneralizedXMLElement("attr", attributes = ["name" => attr[1]],
        content = join(formatted_values, ' '))
end

################################################################################
## Newick
################################################################################
function newickXML(newick::AbstractString;
        using_heights::Bool = true,
        id::String = "startingTree")
    using_dates = !using_heights

    attributes = ["usingHeights" => using_heights, "usingDates" => using_dates]
    return GeneralizedXMLElement("newick", attributes = attributes,
            content = newick, id = id)
end

################################################################################
## treeModel
################################################################################

function treeModelXML(tree::GeneralizedXMLElement, data::DataPairs;
        id::StringOrNothing = nothing)
    children = [tree,
            PassthroughXMLElement("rootHeight",
                parameterXML(id="treeModel.rootHeight")),
            nodeHeightsXML("treeModel.internalNodeHeights",
                internal_nodes = true, root_node = false),
            nodeHeightsXML("treeModel.allInternalNodeHeights",
                internal_nodes = true, root_node = true)]
    for i = 1:length(data)
        trait_name = data[i][1]
        trait_dim = size(data[i][2], 2)
        trait_xml = nodeTraitsXML(trait_name, trait_dim, "$trait_name.leafTraits")

        push!(children, trait_xml)
    end

    attributes = [bn.FIX_HEIGHTS => true]

    return GeneralizedXMLElement("treeModel", children = children, id = id,
                                  attributes = attributes)
end

function nodeHeightsXML(parameter::AbstractGeneralizedXMLElement;
        internal_nodes::Bool = true, root_node::Bool = true)
    attrs = ["internalNodes" => internal_nodes, "rootNode" => root_node]
    return GeneralizedXMLElement("nodeHeights",
        attributes = attrs,
        children=[parameter])
end

function nodeHeightsXML(parameter_id::String; kw_args...)
    return nodeHeightsXML(parameterXML(id=parameter_id); kw_args...)
end

function nodeTraitsXML(trait_name::String, dim::Int,
        parameter::AbstractGeneralizedXMLElement;
        internal_nodes::Bool = false, root_node::Bool = false,
        leaf_nodes::Bool = true,
        as_matrix::Bool = true)
    attrs = ["name" => trait_name,
             "traitDimension" => dim,
             "internalNodes" => internal_nodes,
             "rootNode" => root_node,
             "leafNodes" => leaf_nodes,
             "asMatrix" => as_matrix]
    return GeneralizedXMLElement("nodeTraits",
        attributes = attrs,
        children=[parameter])
end

function nodeTraitsXML(trait_name::String, dim::Int,
        parameter_id::AbstractString; kw_args...)
    return nodeTraitsXML(trait_name, dim, parameterXML(id=parameter_id);
            kw_args...)
end

function find_trait_parameter(tm::GeneralizedXMLElement, trait::String)
    @assert tm.name == "treeModel"
    node_traits = find_element(tm, attributes = Dict("name" => trait))

    return find_element(node_traits, name = "parameter")
end

function transformedTreeModelXML(tm::GeneralizedXMLElement,
                                 param::GeneralizedXMLElement;
                                 id::Nullable{String} = nothing)
    return GeneralizedXMLElement("transformedTreeModel",
        children = [tm, param],
        id = id)
end


################################################################################
## multivariateDiffusionModel
################################################################################

function mbdXML(precision_parameter::AbstractGeneralizedXMLElement;
        id::String = "diffusionModel")
    return GeneralizedXMLElement("multivariateDiffusionModel",
            child=PassthroughXMLElement("precisionMatrix", precision_parameter),
            id = id)
end


################################################################################
## parameters
################################################################################

const ArrayOrNothing{T} = Union{AbstractArray{T}, Nothing} where T
const PARAMETER_NAMES = ("value", "upper", "lower", "dimension")


function parameterXML(;id::StringOrNothing = nothing,
        value::ArrayOrNothing{<:Real} = nothing,
        upper::ArrayOrNothing{<:Real} = nothing,
        lower::ArrayOrNothing{<:Real} = nothing,
        dimension::Union{Int, Nothing} = nothing)

    attributes = Pair{String, Any}[]
    values = (value, upper, lower, dimension)
    for i = 1:length(values)
        add_parameter_attribute!(attributes, PARAMETER_NAMES[i], values[i])
    end

    return GeneralizedXMLElement("parameter", id=id, attributes = attributes)
end

function prior_id(xml::GeneralizedXMLElement)
    return get_id(xml) * ".prior"
end

function add_parameter_attribute!(::Vector{<:Pair{String, <:Any}},
        ::String, ::Nothing)
    # do nothing
end

function add_parameter_attribute!(attributes::Vector{<:Pair{String, <:Any}},
        attr_name::String, attr_value::AbstractArray{<:Real})
    push!(attributes, attr_name => attr_value)
end

function cachedMatrixInverseXML(parameter::AbstractGeneralizedXMLElement;
        id::String = nothing)
    return GeneralizedXMLElement("cachedMatrixInverse", id=id, child=parameter)
end

function compoundSymmetricMatrixXML(
        diagonal_parameter::AbstractGeneralizedXMLElement,
        offdiagonal_parameter::AbstractGeneralizedXMLElement;
        as_correlation::Bool = true,
        is_cholesky::Bool = true,
        strictly_upper::Bool = true,
        id::StringOrNothing = nothing)

    attrs = ["asCorrelation" => as_correlation,
             "isCholesky" => is_cholesky,
             "isStrictlyUpperTriangular" => strictly_upper]
    children = [PassthroughXMLElement("diagonal", diagonal_parameter),
                PassthroughXMLElement("offDiagonal", offdiagonal_parameter)]
    return GeneralizedXMLElement("compoundSymmetricMatrix", attributes = attrs,
            children = children, id = id)
end

function compoundSymmetricMatrixXML(X::AbstractArray{Float64, 2},
        diag_id::String, offdiag_id::String; is_posdef::Bool = true, kw_args...)
    @assert issymmetric(X)
    n = size(X, 1)

    lower_bounds = is_posdef ? zeros(n) : nothing
    diag_parameter = parameterXML(id = diag_id, value = diag(X),
            lower = lower_bounds)

    n_offdiag = div(n * (n - 1), 2)
    offdiag_elements = zeros(n_offdiag)
    ind = 0
    for i = 1:n
        for j = (i + 1):n
            ind += 1
            offdiag_elements[ind] = X[i, j]
        end
    end

    offdiag_parameter = parameterXML(id = offdiag_id, value = offdiag_elements)

    return compoundSymmetricMatrixXML(diag_parameter, offdiag_parameter;
            kw_args...)
end

function compoundSymmetricMatrixXML(dim::Int, diag_id::String,
        offdiag_id::String; kw_args...)
    X = Diagonal(ones(dim))
    return compoundSymmetricMatrixXML(X, diag_id, offdiag_id; kw_args...)
end

function matrixParameterXML(parameters::Vector{<:AbstractGeneralizedXMLElement};
        id::StringOrNothing = nothing)
    return GeneralizedXMLElement("matrixParameter", id = id, children = parameters)
end

function matrixParameterXML(X::AbstractMatrix{Float64};
        id::StringOrNothing = nothing)

    n = size(X, 1)
    ids = isnothing(id) ? fill(nothing, n) : ["$id.$i" for i = 1:n]
    parameters = [parameterXML(id=ids[i], value=X[i, :]) for i = 1:n]
    return matrixParameterXML(parameters, id = id)
end


function maskedParameterXML(parameter::GeneralizedXMLElement,
        mask::AbstractVector{<:Real}; id::Nullable{String} = nothing)
    return GeneralizedXMLElement("maskedParameter", id = id,
            children = [
                    parameter,
                    PassthroughXMLElement("mask", parameterXML(value = mask))
            ])
end


function transformedParameterXML(parameter::GeneralizedXMLElement,
        transform::AbstractGeneralizedXMLElement;
        as_matrix::Bool = false,
        is_multivariate::Bool = as_matrix,
        is_inverse::Bool = false,
        id::Nullable{String} = nothing)
    nm = is_multivariate ? bn.TRANSFORMED_MULTIVARIATE_PARAMETER : bn.TRANSFORMED_PARAMETER
    attrs = [bn.INVERSE => is_inverse]
    if is_multivariate && as_matrix
        attrs = [attrs; bn.AS_MATRIX => as_matrix]
    end
    return GeneralizedXMLElement(nm, children = [parameter, transform],
            attributes = attrs,
            id = id)
end



function compoundParameterXML(params::Vector{<:GeneralizedXMLElement};
                              id::StringOrNothing = nothing)
    return GeneralizedXMLElement(bn.COMPOUND_PARAMETER, id = id, children = params)
end
################################################################################
## traitDataLikelihood
################################################################################
mutable struct TraitLikelihoodOptions
    allowIdentical::Bool
    standardize::Bool
    cacheBranches::Bool
    useTreeLength::Bool
    scaleByTime::Bool
    reportAsMultivariate::Bool
    allowSingular::Bool
end

function TraitLikelihoodOptions(;
        allowIdentical::Bool = true,
        standardize::Bool = true,
        cacheBranches::Bool = true,
        useTreeLength::Bool = false,
        scaleByTime::Bool = true,
        reportAsMultivariate::Bool = true,
        allowSingular::Bool = true)
    return TraitLikelihoodOptions(allowIdentical, standardize, cacheBranches,
            useTreeLength, scaleByTime, reportAsMultivariate, allowSingular)
end

function as_attributes(options::TraitLikelihoodOptions)
    fns = fieldnames(typeof(options))
    n = length(fns)
    attrs = Vector{Pair{String, Any}}(undef, n)
    for i = 1:n
        attrs[i] = string(fns[i]) => getfield(options, fns[i])
    end
    return attrs
end


function get_trait_name(xml::GeneralizedXMLElement)
    if xml.name == "jointPartialsProvider"
        children = [get_trait_name(child) for child in xml.children]
        return join(children, '.') * ".joint"
    end

    return get_attribute(xml, "traitName")
end

function traitDataLikelihoodXML(;
        diffusion_model::GeneralizedXMLElement,
        tree_model::GeneralizedXMLElement,
        extension_model::GeneralizedXMLElement,
        root_mean::Vector{Float64} = zeros(dimension(diffusion_model)),
        prior_sample_size::Float64 = 0.001,
        options = TraitLikelihoodOptions()
        )

    option_attrs = as_attributes(options)
    trait_name = get_trait_name(extension_model)
    push!(option_attrs, "traitName" => trait_name)

    mean_param = PassthroughXMLElement("meanParameter",
            parameterXML(value=root_mean))
    pss_param = PassthroughXMLElement("priorSampleSize",
            parameterXML(value=[prior_sample_size]))

    root_prior = PassthroughXMLElement("conjugateRootPrior",
            [mean_param, pss_param])
    return GeneralizedXMLElement("traitDataLikelihood", id="$trait_name.likelihood",
            attributes = option_attrs,
            children = [tree_model, diffusion_model, extension_model,
                    root_prior]
            )
end

const PARTIALS_PROVIDERS = ["integratedFactorModel",
                            "repeatedMeasuresModel",
                            "jointPartialsProvider"]

function find_partials_provider(xml::GeneralizedXMLElement)
    provider = find_element(xml, names = PARTIALS_PROVIDERS)
end

################################################################################
## Variance priors
################################################################################

function lkjPriorXML(parameter::GeneralizedXMLElement, dimension::Int;
        shape_parameter::Float64 = 1.0,
        id::String = prior_id(parameter)
        )

    attrs = Pair{String, Any}["shapeParameter" => shape_parameter, "dimension" => dimension]
    return GeneralizedXMLElement("LKJCorrelationPrior", id = id,
            attributes = attrs,
            child = PassthroughXMLElement("data", parameter))
end

function halfTPriorXML(parameter::GeneralizedXMLElement;
        df::Int = 1, scale::Float64 = 2.5, id::String = prior_id(parameter))

    return GeneralizedXMLElement("halfTPrior", id = id, child = parameter,
            attributes = Pair{String, Any}["scale" => scale, "df" => df]
    )
end

function determinantPriorXML(parameter::GeneralizedXMLElement;
        shape_parameter::Real = 1.0,
        id::String = prior_id(parameter))

    attrs = Pair{String, Any}["shapeParameter" => shape_parameter]
    return GeneralizedXMLElement("determinantPrior", id = id, attributes = attrs,
        child = parameter)
end


################################################################################
## gradients
################################################################################

function compoundGradientXML(
        sub_gradients::Vector{<:Vector{<:GeneralizedXMLElement}};
        id::StringOrNothing = nothing,
        gradient_ids::Vector{<:StringOrNothing} = fill(nothing, length(sub_gradients))
        )
    n = length(sub_gradients)
    @assert length(gradient_ids) == n
    gradients = [gradientXML(sub_gradients[i]) for i = 1:n]
    return GeneralizedXMLElement("compoundGradient", id = id,
            children = gradients)
end

function gradientXML(children::Vector{<:GeneralizedXMLElement})
    return GeneralizedXMLElement("gradient", children = children)
end

function jointGradientXML(gradients::Vector{<:GeneralizedXMLElement};
        id::StringOrNothing = nothing)
    return GeneralizedXMLElement("jointGradient", id = id, children = gradients)
end

function diffusionGradientXML(;trait_likelihood::GeneralizedXMLElement,
        precision_parameter::GeneralizedXMLElement,
        id::StringOrNothing = nothing,
        parameter::String = "both")

    trait_name = get_attribute(trait_likelihood, "traitName")
    precision_gradient = GeneralizedXMLElement("precisionGradient",
            attributes = ["parameter" => parameter,"traitName" => trait_name],
            children = [trait_likelihood, precision_parameter])
    return GeneralizedXMLElement("diffusionGradient", id = id,
            child = precision_gradient)
end

################################################################################
## hmc
################################################################################

function hmcXML(;
        gradient::GeneralizedXMLElement,
        parameter::GeneralizedXMLElement,
        transform::Union{GeneralizedXMLElement, Nothing} = nothing,
        mask_parameter::Nullable{GeneralizedXMLElement} = nothing,
        mask::Nullable{Vector{<:Real}} = nothing,
        weight::Float64 = 1.0,
        n_steps::Int = 10,
        step_size::Float64 = 0.01,
        draw_variance::Float64 = 1.0,
        gradient_check_count::Int = 0,
        gradient_check_tolerance::Float64 = 0.001,
        is_geodesic::Bool = false,
        orthogonality_structure::Vector{Vector{Int}} = Vector{Int}[])

    if !isnothing(mask_parameter) && !isnothing(mask)
        throw(ArgumentError("cannot supply both 'mask' and 'mask_parameter' keyword arguments"))
    end

    attrs = Pair{String, Any}["weight" => weight, "nSteps" => n_steps,
            "stepSize" => step_size,
            "drawVariance" => draw_variance,
            "gradientCheckCount" => gradient_check_count,
            "gradientCheckTolerance" => gradient_check_tolerance]

    children = AbstractGeneralizedXMLElement[gradient, parameter]
    if !isnothing(transform)
        push!(children, transform)
    end

    if !isnothing(mask)
        mask_parameter = parameterXML(value = mask)
    end

    if !isnothing(mask_parameter)
        push!(children, PassthroughXMLElement("mask", mask_parameter))
    end

    hmc_name = is_geodesic ? "geodesicHamiltonianMonteCarloOperator" :
            "hamiltonianMonteCarloOperator"

    if is_geodesic && length(orthogonality_structure) > 0
        structures = [GeneralizedXMLElement("x", attributes = ["rows" => rows])
                for rows in orthogonality_structure]
        push!(children, PassthroughXMLElement("orthogonalityStructure", structures))
    end

    return GeneralizedXMLElement(hmc_name, children = children,
            attributes = attrs)
end

function lkjTransformXML(dim::Int)
    return GeneralizedXMLElement("LKJTransform",
            attributes = ["dimension" => dim])
end

function transformXML(type::String; dim::Nullable{Int} = nothing)
    return GeneralizedXMLElement("transform",
            attributes = ["type" => type, "dim" => dim])
end

function multivariateCompoundTransformXML(transforms::GeneralizedXMLElement...)
    return GeneralizedXMLElement("multivariateCompoundTransform",
            children = collect(transforms))
end


################################################################################
## other operators
################################################################################
"""
<normalGammaPrecisionGibbsOperator weight="1.0">
    <prior>
        <gammaPrior idref="factorPrecision.prior"/>
    </prior>
    <normalExtension treeTraitName="factors">
        <integratedFactorModel idref="factorModel"/>
        <traitDataLikelihood idref="traitLikelihood"/>
    </normalExtension>
</normalGammaPrecisionGibbsOperator>
"""
function normalGammaGibbsXML(;prior::GeneralizedXMLElement,
        provider::GeneralizedXMLElement,
        weight::Float64 = 1.0)
    return GeneralizedXMLElement("normalGammaPrecisionGibbsOperator",
            children = [PassthroughXMLElement("prior", prior), provider],
            attributes = ["weight" => weight])
end

function normalExtensionXML(;
        extension::GeneralizedXMLElement,
        likelihood::GeneralizedXMLElement,
        trait_name::String = get_attribute(likelihood, "traitName"))
    return GeneralizedXMLElement("normalExtension",
            children = [extension, likelihood],
            attributes = ["treeTraitName" => trait_name])
end

"""
<randomWalkOperator windowSize="1.0" weight="3">
    <parameter idref="exponential.growthRate"/>
</randomWalkOperator>
"""
function randomWalkXML(parameter::GeneralizedXMLElement;
        window_size::Float64 = 1.0, weight::Float64 = 1.0)
    return GeneralizedXMLElement("randomWalkOperator", child = parameter,
            attributes = [bn.WINDOW_SIZE => window_size, bn.WEIGHT => weight])
end

"""
<transformedParameterRandomWalkOperator windowSize="0.1" weight="10" checkValid="true">
    <transformedMultivariateParameter idref="corr"/>
    <correlationBounds dimension="4"/>
</transformedParameterRandomWalkOperator>
"""

function transformedRandomWalkXML(parameter::GeneralizedXMLElement;
        bounds::Nullable{<:GeneralizedXMLElement} = nothing,
        weight::Real = 1,
        check_valid::Bool = true,
        window_size::Real = 0.1,
        mask::Vector{Int} = Int[]
        )
    attrs = [bn.WEIGHT => weight, bn.WINDOW_SIZE => window_size, bn.CHECK_VALID => check_valid]
    children = AbstractGeneralizedXMLElement[parameter]
    if !isnothing(bounds)
        push!(children, bounds)
    end

    if length(mask) > 0
        mask_xml = PassthroughXMLElement(bn.UPDATE_INDEX, parameterXML(value = mask))
        push!(children, mask_xml)
    end

    return GeneralizedXMLElement(bn.TRANSFORMED_RANDOM_WALK,
            children = children,
            attributes= attrs)
end



################################################################################
## distributions
################################################################################

"""
<distributionLikelihood id="L.prior">
    <data>
        <matrixParameter idref="L"/>
    </data>
    <distribution>
        <normalDistributionModel>
            <mean>
                <parameter value="0.0"/>
            </mean>
            <stdev>
                <parameter value="1.0" lower="0"/>
            </stdev>
        </normalDistributionModel>
    </distribution>
</distributionLikelihood>
"""

function distributionLikelihoodXML(;data_parameter::GeneralizedXMLElement,
        distribution_model::GeneralizedXMLElement,
        id::String = prior_id(data_parameter))
    return GeneralizedXMLElement("distributionLikelihood", id = id,
            children = [PassthroughXMLElement("data", data_parameter),
                    PassthroughXMLElement("distribution", distribution_model)])
end

function normalDistributionModelXML(;mean_parameter::GeneralizedXMLElement,
        stdev_parameter::GeneralizedXMLElement)
    return GeneralizedXMLElement("normalDistributionModel",
            children = [PassthroughXMLElement("mean", mean_parameter),
                    PassthroughXMLElement("stdev", stdev_parameter)])
end

function standardNormalDistributionXML()
    return normalDistributionModelXML(mean_parameter = parameterXML(value=[0]),
            stdev_parameter = parameterXML(value=[1], lower=[0]))
end

"""
<gammaPrior id="factorPrecision.prior" scale="1.0" shape="1.0">
    <parameter idref="factorPrecision"/>
</gammaPrior>
"""

function gammaPriorXML(parameter::GeneralizedXMLElement;
        scale::Float64 = 1.0, shape::Float64 = 1.0,
        id::String = prior_id(parameter))
    return GeneralizedXMLElement("gammaPrior", id = id, child = parameter,
            attributes = ["scale" => scale, "shape" => shape])
end

################################################################################
## mcmc
################################################################################

function mcmcXML(joint::GeneralizedXMLElement, operators::GeneralizedXMLElement,
        screen_log::GeneralizedXMLElement,
        file_log::GeneralizedXMLElement,
        chain_length::Int;
        auto_optimize::Bool = true)
    return GeneralizedXMLElement("mcmc", id = "mcmc",
        children = [joint, operators, screen_log, file_log],
        attributes = Pair{String, Any}["chainLength" => chain_length,
                "autoOptimize" => auto_optimize])
end

function mcmcXML(org::Organizer, operators::GeneralizedXMLElement;
            chain_length::Int,
            screen_logEvery::Int,
            file_logEvery::Int,
            file_name::String,
            overwrite::Bool = true)
    likelihood = GeneralizedXMLElement("likelihood", id = "likelihood",
            children = org.likelihoods)
    prior = GeneralizedXMLElement("prior", id = "prior", children = org.priors)
    joint = GeneralizedXMLElement("joint", id="joint",
            children = [likelihood, prior])

    default_loggables = [joint, likelihood, prior]
    screen_log = logXML(default_loggables, screen_logEvery)
    file_log = logXML([default_loggables; org.loggables], file_logEvery,
        file_name = file_name,
        overwrite = overwrite)

    return mcmcXML(joint, operators, screen_log, file_log, chain_length)
end

function mcmcXML(org::Organizer, operators::GeneralizedXMLElement,
        mcmc_options::MCMCOptions; file_name::String)

    return mcmcXML(org, operators, chain_length = mcmc_options.chain_length,
        file_logEvery = mcmc_options.file_log_every,
        screen_logEvery = mcmc_options.screen_log_every, file_name = file_name)
end

function logXML(elements::Vector{<:GeneralizedXMLElement}, log_every::Int;
        id::StringOrNothing = nothing,
        overwrite::Union{Nothing, Bool} = nothing,
        file_name::StringOrNothing = nothing)

    attrs = ["fileName" => file_name,
            "logEvery" => log_every,
            "overwrite" => overwrite]
    return GeneralizedXMLElement("log", id = id, children = elements,
            attributes = attrs)
end


################################################################################
## Latent Liability
################################################################################

"""
<generalDataType id="discreteStates">
    <state code="0"/>
    <state code="1"/>
    <ambiguity code="?" states="01"/>
</generalDataType>
"""

function generalDataTypeXML(;
        id::String,
        states::Vector,
        ambiguity_code::String = "?")
    states_xml = [GeneralizedXMLElement(bn.STATE, attributes = [bn.CODE => state])
                for state in states]
    push!(states_xml, GeneralizedXMLElement(bn.AMBIGUITY, attributes = [
                                                    bn.CODE => ambiguity_code,
                                                    bn.STATES => join(states)]))
    return GeneralizedXMLElement(bn.GENERAL_DATA_TYPE, children = states_xml, id = id)
end


"""
<alignment id="traits.alignment">
    <generalDataType idref="discreteStates"/>
    <sequence><taxon idref="139_TheraceraPennigera_1_NF_N_RAW_AOB"/>1
    </sequence>
    <sequence><taxon idref="127_TayuvaLilacina_8_NP_N_RAW_AOB"/>0
    </sequence>
    <sequence><taxon idref="127_TayuvaLilacina_7_NP_N_RAW_AOB"/>0
    </sequence>
</alignment>
"""

function alignmentXML(taxaElements::Vector{GeneralizedXMLElement},
                      dataType::GeneralizedXMLElement,
                      sequences::Matrix,
                      taxa::Vector{<:AbstractString};
                      id::StringOrNothing = nothing)

    n = length(taxa)
    @assert size(sequences, 1) == length(taxaElements) == n
    seqs = Vector{GeneralizedXMLElement}(undef, n)
    for i = 1:n
        taxon = taxa[i]
        @assert get_id(taxaElements[i]) == taxon
        seqXML = GeneralizedXMLElement(bn.SEQUENCE,
                    children = [taxaElements[i]],
                    content = join(Int.(sequences[i, :]), ' '))
        seqs[i] = seqXML
    end

    return GeneralizedXMLElement(bn.ALIGNMENT, children = [dataType; seqs], id = id)
end


"""
<patterns id="traits.patterns" from="1" unique="false">
    <alignment idref="traits.alignment"/>
</patterns>
"""

function patternsXML(alignment::GeneralizedXMLElement;
                     id::StringOrNothing = get_id(alignment) * ".patterns")
    attrs = [bn.FROM => "1", bn.UNIQUE => false]
    return GeneralizedXMLElement(bn.PATTERNS, children = [alignment],
                                  attributes = attrs, id = id)
end



"""
<orderedLatentLiabilityLikelihood NTraits="1" NData="1" id="traits.latentLiability">
    <patterns idref="traits.patterns"/>
    <treeModel idref="treeModel"/>
    <tipTrait>
        <parameter idref="leafTraits"/>
    </tipTrait>
    <threshold>
        <compoundParameter id="leafTraits.threshold">
            <parameter value="0.0"/>
        </compoundParameter>
    </threshold>
    <numClasses>
        <parameter value="2.0"/>
    </numClasses>
</orderedLatentLiabilityLikelihood>
"""

function orderedLatentLiabilityLikelihoodXML(;
            patterns::GeneralizedXMLElement,
            treeModel::GeneralizedXMLElement,
            trait_parameter::GeneralizedXMLElement,
            id::StringOrNothing = nothing
            )
    @warn "latent liability currently only implemented for a 1-D trait"
    tip_child = PassthroughXMLElement(bn.TIP_TRAIT, trait_parameter)
    thresh_param = compoundParameterXML([parameterXML(value = [0.0])])
    threshold = PassthroughXMLElement(bn.THRESHOLD, thresh_param)
    num_classes = PassthroughXMLElement(bn.NUM_CLASSES, parameterXML(value=[2.0]))

    children = [patterns, treeModel, tip_child, threshold, num_classes]
    attrs = [bn.N_TRAITS => "1", bn.N_DATA => "1"]
    return GeneralizedXMLElement(bn.ORDERED_LATENT_LIABILITY,
            children = children, attributes = attrs, id = id)
end


"""
<extendedLatentLiabilityGibbsOperator weight="1.0" traitName="">
    <traitDataLikelihood idref="traitLikelihood"/>
    <orderedLatentLiabilityLikelihood idref="traits.latentLiability"/>
</extendedLatentLiabilityGibbsOperator>
"""

function extendedLatentLiabilityGibbsOperatorXML(;
        likelihood::GeneralizedXMLElement,
        latent_liability::GeneralizedXMLElement,
        trait_name::AbstractString,
        weight::Real = 1
        )
    return GeneralizedXMLElement(bn.EXTENDED_GIBBS_OPERATOR,
            children = [likelihood, latent_liability],
            attributes = [bn.WEIGHT => weight, bn.TRAIT_NAME => trait_name])
end


################################################################################
## Statistics
################################################################################


"""
<traitLogger id="traitLogger" traitName="factors" taxonNameExplicit="true" nodes="external">
    <treeModel idref="treeModel"/>
    <traitDataLikelihood idref="traitLikelihood"/>
</traitLogger>
"""

function traitLoggerXML(;
        tree_model::GeneralizedXMLElement,
        trait_likelihood::GeneralizedXMLElement,
        trait_name::String = get_trait_name(trait_likelihood),
        id::String = trait_name * ".traitLogger",
        taxon_name_explicit::Bool = true,
        nodes::String = bn.EXTERNAL)
    return GeneralizedXMLElement(bn.TRAIT_LOGGER, id = id,
            children = [tree_model, trait_likelihood],
            attributes = [bn.TRAIT_NAME => trait_name,
                    bn.TAXON_EXPLICIT => taxon_name_explicit,
                    bn.NODES => nodes])
end


"""
<factorProportionStatistic id="factorProportion">
    <integratedFactorModel idref="factorModel"/>
    <traitDataLikelihood idref="traitLikelihood"/>
</factorProportionStatistic>
"""

function factorProportionStatisticXML(;
        factor_model::GeneralizedXMLElement,
        trait_likelihood::GeneralizedXMLElement,
        id::String = get_attribute(factor_model, bn.TRAIT_NAME) * ".proportion")
    return GeneralizedXMLElement(bn.FACTOR_PROPORTION_STATISTIC, id = id,
            children = [factor_model, trait_likelihood])
end


"""
<blombergsK id="kstat" traitName="X">
    <traitDataLikelihood idref="traitLikelihood"/>
</blombergsK>
"""

function blombergsKStatisticXML(;
        trait_likelihood::GeneralizedXMLElement,
        trait_name::AbstractString = get_trait_name(trait_likelihood),
        id::AbstractString = "K.$trait_name")

    return GeneralizedXMLElement("blombergsK",
        children = [trait_likelihood],
        attributes = [bn.TRAIT_NAME => trait_name],
        id = id)
end



function postOrderRootMeanXML(;
        trait_likelihood::GeneralizedXMLElement,
        id::AbstractString = "rootMean.$(get_trait_name(trait_likelihood))")
    return GeneralizedXMLElement("postOrderRootMeanStatistic",
            children = [trait_likelihood],
            id = id)
end

################################################################################
## Adjacency
################################################################################


"""
<adjacentTaxa>
    <treeModel idref="treeModel"/>
    <adjacency sharedPerimeter="1">
        <taxonData perimeter="4" area="2">
            <taxon idref="02013"/>
        </taxonData>
        <taxonData perimeter="4" area="2">
            <taxon idref="02016"/>
        </taxonData>
    </adjacency>
    ...
</adjacentTaxa>

"""

function taxonDataXML(taxon::GeneralizedXMLElement;
                      perimeter::Float64,
                      area::Float64)
    return GeneralizedXMLElement("taxonData",
                                 children = [taxon],
                                 attributes = ["perimeter" => perimeter,
                                               "area" => area])
end

function adjacentTaxaXML(;
                         taxa::Vector{GeneralizedXMLElement},
                         tree_model::GeneralizedXMLElement,
                         adjacent_taxa::DataFrame,
                         id::AbstractString = "adjacentTaxa",
                         )
    tdos = [
            taxonDataXML(
                    find_taxon(taxa,
                               adjacent_taxa.orig[i]),
                               perimeter = adjacent_taxa.orig_perim[i],
                               area = adjacent_taxa.orig_area[i])
            for i = 1:nrow(adjacent_taxa)]
    tdds = [
            taxonDataXML(
                    find_taxon(taxa,
                                adjacent_taxa.dest[i]),
                                perimeter = adjacent_taxa.dest_perim[i],
                                area = adjacent_taxa.dest_area[i])
            for i = 1:nrow(adjacent_taxa)]
    adjacencies = [PassthroughXMLElement("adjacency",
            [tdos[i], tdds[i]])
            for i = 1:nrow(adjacent_taxa)]

    return GeneralizedXMLElement("adjacentTaxa", children = [tree_model; adjacencies], id = id)
end

function find_taxon(taxa::Vector{GeneralizedXMLElement}, id::AbstractString)
    taxa[findfirst(x -> x.id == id, taxa)]
end