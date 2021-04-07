# getting elements
function get_id(el::MyXMLElement)
    return el.id
end

function get_data(bx::BEASTXMLElement)
    return find_element(bx, DataXMLElement)
end

function get_newick(bx::BEASTXMLElement)
    return find_element(bx, NewickXMLElement)
end

function get_treeModel(bx::BEASTXMLElement)
    return find_element(bx, TreeModelXMLElement)
end

function get_mbd(bx::BEASTXMLElement)
    return find_element(bx, MBDXMLElement)
end

function get_extension(bx::BEASTXMLElement)
    return find_element(bx, ModelExtensionXMLElement)
end

function get_repeatedMeasures(bx::BEASTXMLElement)
    return find_element(bx, RepeatedMeasuresXMLElement)
end

function get_integratedFactorModel(bx::BEASTXMLElement)
    return find_element(bx, IntegratedFactorsXMLElement)
end

function get_latentFactorModel(bx::BEASTXMLElement)
    return find_element(bx, LatentFactorModelXMLElement)
end

function get_traitLikelihood(bx::BEASTXMLElement)
    return find_element(bx, TraitLikelihoodXMLElement)
end

function get_operators(bx::BEASTXMLElement)
    ops = find_element(bx, OperatorsXMLElement)
    return ops.els
end

function get_mcmc(bx::BEASTXMLElement)
    return find_element(bx, MCMCXMLElement)
end

function get_timer(bx::BEASTXMLElement)
    return find_element(bx, TimerXMLElement)
end

function get_loadings_op(bx::BEASTXMLElement; component::String = "")
    loadings = get_loadings(bx)
    if typeof(loadings) <: ScaledOrthogonalMatrix
        if component == "scale"
            loadings = loadings.scale
        elseif component == "matrix"
            loadings = loadings.U
        else
            throw(ArgumentError("The loadings is a " * name(loadings) *
                " with components \"scale\" and \"matrix\"." *
                " You have set component=\"$component\"."))
        end
    end

    ops = get_operators(bx)
    for op in ops
        if get_parameter(op) === loadings
            return op
        end
    end
    error("No loadings operator.")
end

function get_op_by_parameter(bx::BEASTXMLElement, parameter::MyXMLElement)
    ops = get_operators(bx)
    for op in ops
        if get_parameter(op) === parameter
            return op
        end
    end
    error("cannot find operator for parameter with id '$(get_id(parameter))'")
end

function get_multiplicative_gamma_op(bx::BEASTXMLElement)
    ops = get_operators(bx)
    for op in ops
        t = typeof(op)
        if t <: NormalGammaPrecisionOperatorXMLElement
            if typeof(op.ggp) <: MultipilcativeGammaGibbsProvider
                return op
            end
        end
    end
    error("No multiplicative gamma gibbs operator.")
end

function get_loadings_scale(bx::BEASTXMLElement)
    ifm = get_integratedFactorModel(bx)
    return get_loadings_scale(ifm)
end

function get_factor_model(bx::BEASTXMLElement)
    fac_models = [IntegratedFactorsXMLElement, LatentFactorModelXMLElement]
    for model in fac_models
        fac_model = find_element(bx, model)
        if !isnothing(fac_model)
            return fac_model
        end
    end

    return nothing
end

function get_loadings(bx::BEASTXMLElement)
    fac_model = get_factor_model(bx)
    loadings = get_loadings(fac_model)
    return loadings
end

function get_factor_precisions(bx::BEASTXMLElement)
    fac_model = get_factor_model(bx)
    precs = get_precision(fac_model)
    return precs
end