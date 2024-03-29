
mutable struct CompoundLikelihood <: MyXMLElement
    el::XMLOrNothing
    id::String
    likelihoods::Vector{<:MyXMLElement}
    already_made::Vector{Bool}
    function CompoundLikelihood(id::String, likelihoods::Vector{<:MyXMLElement};
                                already_made::Vector{Bool} = fill(true, length(likelihoods)))

        @assert length(already_made) == length(likelihoods)
        return new(nothing, id, likelihoods, already_made)
    end
end

function name(cl::CompoundLikelihood)
    return bn.LIKELIHOOD
end

function get_loggables(cl::CompoundLikelihood)
    make_xml(cl)
    return cl.el
end

function add_likelihood!(cp::CompoundLikelihood, xml::MyXMLElement; already_made::Bool = true)
    push!(cp.likelihoods, xml)
    push!(cp.already_made, already_made)
end

mutable struct MCMCXMLElement <: MyXMLElement
    el::XMLOrNothing
    likelihoods::CompoundLikelihood
    priors::Vector{<:MyXMLElement}
    operators::OperatorsXMLElement
    chain_length::Int
    file_logEvery::Int
    screen_logEvery::Int
    loggables::LoggablesXMLElement
    filename::String
    log_files::Bool
    attrs::Dict{String, String}
    extras::Vector{<:MyXMLElement}

    function MCMCXMLElement(tl_el::TraitLikelihoodXMLElement,
                    mbd_el::MBDXMLElement,
                    rm_el::RepeatedMeasuresXMLElement,
                    os_el::OperatorsXMLElement; chain_length = 100)
        fle = max(1, div(chain_length, 1000))
        sle = 1000
        filename = "defaultFile"
        lg_el = LoggablesXMLElement([mbd_el, rm_el], [true, true])
        attrs = Dict(bn.AUTO_OPTIMIZE => bn.TRUE)

        return new(nothing, CompoundLikelihood("likelihood", [tl_el]), [mbd_el, rm_el], os_el, chain_length,
                    fle, sle, lg_el, filename, true, attrs,
                    MyXMLElement[])
    end

    function MCMCXMLElement(tl_el::TraitLikelihoodXMLElement,
                    mbd_el::MBDXMLElement,
                    os_el::OperatorsXMLElement; chain_length = 100)
        fle = max(1, div(chain_length, 1000))
        sle = 1000
        filename = "defaultFile"
        lg_el = LoggablesXMLElement([mbd_el], [true])
        attrs = Dict(bn.AUTO_OPTIMIZE => bn.TRUE)

        return new(nothing, CompoundLikelihood("likelihood", [tl_el]), [mbd_el], os_el, chain_length,
                    fle, sle, lg_el, filename, true, attrs,
                    MyXMLElement[])
    end

    function MCMCXMLElement(tl_el::TraitLikelihoodXMLElement,
                    mbd_el::MBDXMLElement,
                    rm_el::RepeatedMeasuresXMLElement,
                    os_el::OperatorsXMLElement,
                    loggables::LoggablesXMLElement;
                    chain_length = 100)
        fle = max(1, div(chain_length, 1000))
        sle = 1000
        filename = "defaultFile"
        attrs = Dict(bn.AUTO_OPTIMIZE => bn.TRUE)


        return new(nothing, [tl_el], CompoundLikelihood("likelihood", [mbd_el, rm_el]), os_el, chain_length,
                    fle, sle, loggables, filename, true, attrs,
                    MyXMLElement[])
    end

    function MCMCXMLElement(tl_el::TraitLikelihoodXMLElement,
                    mbd_el::MBDXMLElement,
                    if_el::Union{IntegratedFactorsXMLElement, LatentFactorModelXMLElement},
                    os_el::OperatorsXMLElement; chain_length = 100)
        fle = max(1, div(chain_length, 1000))
        sle = 1000
        filename = "defaultFile"
        lg_el = LoggablesXMLElement([mbd_el, if_el], [true, true])
        attrs = Dict(bn.AUTO_OPTIMIZE => bn.TRUE)

        return new(nothing, CompoundLikelihood("likelihood", [tl_el, if_el]), [mbd_el, if_el], os_el, chain_length,
                    fle, sle, lg_el, filename, true, attrs,
                    MyXMLElement[])
    end
end

function set_options!(mcmc::MCMCXMLElement, options::MCMCOptions)
    mcmc.chain_length = options.chain_length
    mcmc.file_logEvery = options.file_log_every
    mcmc.screen_logEvery = options.screen_log_every
    set_full_eval!(mcmc, options.likelihood_check_count)
end

function set_chain_length!(mcmc::MCMCXMLElement, cl::Int)
    mcmc.chain_length = cl
end

function set_full_eval!(mcmc::MCMCXMLElement, n_eval::Int)
    mcmc.attrs[bn.FULL_EVALUATION] = string(n_eval)
end

function make_xml(mc_el::MCMCXMLElement)
    el = new_element(bn.MCMC)
    set_attribute(el, bn.ID, bn.MCMC)

    mc_el.attrs[bn.CHAINLENGTH] = string(mc_el.chain_length)
    set_attributes(el, mc_el.attrs)

    posterior_el = new_child(el, bn.POSTERIOR)
    set_attribute(posterior_el, bn.ID, bn.POSTERIOR)
    prior_el = new_child(posterior_el, bn.PRIOR)
    set_attribute(prior_el, bn.ID, bn.PRIOR)

    for prior in mc_el.priors
        add_ref_els(prior_el, get_priors(prior))
    end
    like_el = make_xml(mc_el.likelihoods)
    add_child(posterior_el, like_el)

    add_ref_el(el, mc_el.operators.el)

    screen_log_el = new_child(el, bn.LOG)
    attrs = [(bn.ID, bn.SCREEN_LOG_ID), (bn.LOGEVERY, string(mc_el.screen_logEvery))]
    set_attributes(screen_log_el, attrs)

    default_logs = [posterior_el, prior_el, like_el]
    add_screenLog_cols(screen_log_el, default_logs)

    if mc_el.log_files

        file_log_el = new_child(el, bn.LOG)
        attrs = [(bn.ID, bn.FILE_LOG_ID), (bn.LOGEVERY, string(mc_el.file_logEvery)),
                (bn.FILENAME, "$(mc_el.filename).log"), (bn.OVERWRITE, bn.TRUE)]
        set_attributes(file_log_el, attrs)

        file_logs = Vector{XMLElement}(undef, 0)
        for myxml in mc_el.loggables.els
            file_logs = [file_logs; get_loggables(myxml)]
        end
        # file_logs = [get_loggable(loggable_el) for loggable_el in mc_el.loggables.els]

        final_file_logs = [default_logs; file_logs]
        for log in final_file_logs
            add_ref_el(file_log_el, log)
        end

    end

    for ext_el in mc_el.extras
        make_xml(ext_el)
        add_child(el, ext_el.el)
    end

    mc_el.el = el
    return el
end

function add_screenLog_cols(pel::XMLElement, cols::Array{XMLElement};
            width::Int = bn.WIDTH_DEFAULT, dp::Int = bn.DP_DEFAULT)
    for el in cols
        label = attribute(el, bn.ID, required = true)
        col_el = new_element(bn.COLUMN)
        attrs = [(bn.LABEL, label), (bn.DP, string(dp)), (bn.WIDTH, string(width))]
        set_attributes(col_el, attrs)
        add_ref_el(col_el, el)
        add_child(pel, col_el)
    end
end

function set_screen_logEvery(mcmc::MCMCXMLElement, sle::Int)
    mcmc.screen_logEvery = sle
end

function set_file_logEvery(mcmc::MCMCXMLElement, fle::Int)
    mcmc.file_logEvery = fle
end

function set_filename(mcmc::MCMCXMLElement, filename::String)
    mcmc.filename = filename
end

function merge_mcmc!(mcmc::MCMCXMLElement, pmcmc::ParsedMCMCXMLElement;
                     separate_logs::Bool = false)

    trait_id = "traitLikelihoods"
    trait_likelihood = mcmc.likelihoods
    trait_likelihood.id = trait_id


    mcmc.likelihoods = CompoundLikelihood("likelihood",
                            [mcmc.likelihoods; pmcmc.likelihoods],
                            already_made = [false; fill(true, length(pmcmc.likelihoods))])
    add_loggable(mcmc.loggables, trait_likelihood, already_made = true)

    mcmc.priors = [mcmc.priors; pmcmc.priors]

    if separate_logs
        push!(mcmc.extras, pmcmc.file_log)
    else
        push!(mcmc.loggables, LoggablesXMLElement(pmcmc.file_loggables))
    end

    push!(mcmc.extras, pmcmc.tree_log)
end

function add_prior!(mcmc::MCMCXMLElement, xml::MyXMLElement)
    push!(mcmc.priors, xml)
end

function add_likelihood!(mcmc::MCMCXMLElement, xml::MyXMLElement)
    add_likelihood!(mcmc.likelihoods, xml)
end

function add_loggable!(mcmc::MCMCXMLElement, xml::MyXMLElement)
    add_loggable(mcmc.loggables, xml)
end


function make_xml(cl::CompoundLikelihood)
    el = new_element(bn.LIKELIHOOD)
    set_attribute(el, bn.ID, cl.id)
    for i = 1:length(cl.likelihoods)
        if cl.already_made[i]
            add_ref_el(el, cl.likelihoods[i])
        else
            c_el = make_xml(cl.likelihoods[i])
            add_child(el, c_el)
        end
    end
    cl.el = el
    return el
end


