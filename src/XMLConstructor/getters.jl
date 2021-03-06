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

function get_loadings_op(bx::BEASTXMLElement)
    loadings = get_loadings(bx)
    ops = get_operators(bx)
    for op in ops
        if get_parameter(op) === loadings
            return op
        end
    end
    error("No loadings operator.")
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

function get_loadings(bx::BEASTXMLElement)
    fac_models = [IntegratedFactorsXMLElement, LatentFactorModelXMLElement]
    fac_model = nothing
    for model in fac_models
        fac_model = find_element(bx, model)
        if !isnothing(fac_model)
            break
        end
    end
    loadings = get_loadings(fac_model)
    return loadings
end