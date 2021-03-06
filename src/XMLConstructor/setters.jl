## set parts of BeastXML

function set_chain_length!(bx::BEASTXMLElement, cl::Int)
    mcmc_el = get_mcmc(bx)
    set_chain_length!(mcmc_el, cl)
end

function set_screen_logEvery(bx::BEASTXMLElement, sle::Int)
    mcmc_el = get_mcmc(bx)
    set_screen_logEvery(mcmc_el, sle)
end

function set_file_logEvery(bx::BEASTXMLElement, fle::Int)
    mcmc_el = get_mcmc(bx)
    set_file_logEvery(mcmc_el, fle)
end

function set_filename(bx::BEASTXMLElement, filename::AbstractString)
    mcmc_el = get_mcmc(bx)
    if isnothing(mcmc_el)
        error("No MCMC component. Cannot set log file names.")
    end
    set_filename(mcmc_el, string(filename))
    timer_el = get_timer(bx)
    if !isnothing(timer_el)
        set_filename(timer_el, "$(filename)_timer.txt")
    end
end

function set_data_dates(bx::BEASTXMLElement, dates::AbstractVector{Float64})
    data_el = get_data(bx)
    set_dates(data_el, dates)
end

function set_full_eval!(bx::BEASTXMLElement, n_eval::Int)
    mcmc = get_mcmc(bx)
    set_full_eval!(mcmc, n_eval)
end

## Model Parameters
function set_mbd_precision(bx::BEASTXMLElement, mat::AbstractArray{Float64, 2})
    mbd = get_mbd(bx)
    set_precision(mbd, mat)
end

function set_residual_precision(bx::BEASTXMLElement, mat::AbstractArray{Float64, 2})
    ext = get_extension(bx)
    set_precision(ext, mat)
end

function set_loadings(bx::BEASTXMLElement, L::AbstractArray{Float64, 2})
    loadings = get_loadings(bx)
    set_mat!(loadings, L)
end

function set_muliplicative_gamma_indices(bx::BEASTXMLElement, inds::Vector{Int})
    op = get_multiplicative_gamma_op(bx)
    set_indices!(op, inds)
end

function set_options!(bx::BEASTXMLElement, options::MCMCOptions)
    mcmc = get_mcmc(bx)
    set_options!(mcmc, options)
end

