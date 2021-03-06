################################################################################
## Parameter models
################################################################################

abstract type ParamsProvider end
# interface:
#   tip_dimension(provider::ParamsProvider)
#   data_dimension(provider::ParamsProvider)
#   make_xmlelement(provider::ParamsProvider, tm::TreeModelXMLElement; ind::Int = 1)

const DEFAULT_PSS = 0.001

function data_dimension(provider::ParamsProvider)
    return tip_dimension(provider)
end


mutable struct DiffusionModel <: ParamsProvider
    Σ::AbstractArray{Float64, 2} # diffusion variance
    μ::AbstractVector{Float64} # diffusion mean
    pss::Float64 # prior sample size of diffusion mean

    function DiffusionModel(Σ::AbstractMatrix{Float64},
                            μ::AbstractVector{Float64},
                            pss::Float64)
        check_posdef(Σ)
        p = size(Σ, 1)
        if length(μ) != p
            error("The mean and variance have incompatible dimensions.")
        end
        if pss <= 0.0
            error("The prior sample size must be positive.")
        end
        return new(Σ, μ, pss)
    end
end

function DiffusionModel(Σ::AbstractMatrix{Float64})
    p = size(Σ, 1)
    return DiffusionModel(Σ, zeros(p), DEFAULT_PSS)
end

function DiffusionModel(p::Int)
    return DiffusionModel(Diagonal(ones(p)), zeros(p), DEFAULT_PSS)
end

function tip_dimension(provider::DiffusionModel)
    return length(provider.μ)
end

mutable struct IntegratedFactorModel <: ParamsProvider
    L::AbstractMatrix{Float64} # KxP loadings matrix
    λ::AbstractVector{Float64} # P-vector of residual variances

    function IntegratedFactorModel(L::AbstractMatrix{Float64},
                                   λ::AbstractVector{Float64})
        k, p = size(L)
        if length(λ) != p
            error("The loadings and residual variances have incompatible " *
                  "dimensions.")
        end
        return new(L, λ)
    end
end

function IntegratedFactorModel(k::Int, p::Int)
    return IntegratedFactorModel(zeros(k, p), ones(p))
end

function tip_dimension(provider::IntegratedFactorModel)
    return size(provider.L, 1)
end

function data_dimension(provider::IntegratedFactorModel)
    return size(provider.L, 2)
end

mutable struct ResidualVarianceModel <: ParamsProvider
    Γ::AbstractMatrix{Float64} # residual variance

    function ResidualVarianceModel(Γ::AbstractMatrix{Float64})
        check_posdef(Γ)
        return new(Γ)
    end
end

function ResidualVarianceModel(p::Int)
    return ResidualVarianceModel(Diagonal(ones(p)))
end

function tip_dimension(provider::ResidualVarianceModel)
    return size(provider.Γ, 1)
end


mutable struct JointProcessModel
    diffusion_model::DiffusionModel
    extensions::AbstractArray{<:ParamsProvider}

    function JointProcessModel(diff::DiffusionModel,
                               extensions::AbstractArray{<:ParamsProvider})
        p_diff = tip_dimension(diff)
        p_ext = 0
        for ext in extensions
            p_ext += tip_dimension(ext)
        end
        if p_diff != p_ext
            error("Incompatible dimensions between diffusion process and " *
                  "extensions")
        end
        return new(diff, extensions)
    end
end

function JointProcessModel(extensions::AbstractArray{<:ParamsProvider})
    p_ext = 0
    for ext in extensions
        p_ext += tip_dimension(ext)
    end
    return JointProcessModel(DiffusionModel(p_ext), extensions)
end

function tip_dimension(jpm::JointProcessModel)
    return tip_dimension(jpm.diffusion_model)
end

function data_dimension(jpm::JointProcessModel)
    dim = 0
    for ext in jpm.extensions
        dim += data_dimension(ext)
    end
    return dim
end



################################################################################
## Data models
################################################################################

const DEFAULT_TRAIT_NAME = "traits"

mutable struct DataModel
    taxa::AbstractVector{<:AbstractString}
    data::Matrix{Float64}
    partitions::Vector{UnitRange{Int}}
    trait_names::Vector{<:AbstractString}

    function DataModel(taxa::AbstractVector{<:AbstractString},
                       data::AbstractMatrix{Float64}, name::AbstractString)
        p = size(data, 2)
        return new(taxa, Matrix(data), [1:p], [name])
    end

    function DataModel(taxa::AbstractVector{<:AbstractString},
                       data::AbstractArray{<:AbstractArray{Float64}},
                       trait_names::AbstractArray{<:AbstractString})
        if length(data) != length(trait_names)
            error("The number of data matrices must equal the number of " *
                  "trait names")
        end
        n_mats = length(data)
        n = length(taxa)
        partitions = Vector{UnitRange}(undef, n_mats)
        stop = 0
        for i = 1:n_mats
            start = stop + 1
            n_sub, p = size(data[i])
            stop += p
            if n_sub != n
                error("All matrices must have the same number of rows as taxa.")
            end
            partitions[i] = start:stop
        end
        p_total = stop
        mat = zeros(n, p_total)
        for i = 1:n_mats
            mat[1:n, partitions[i]] .= data[i]
        end

        return new(taxa, mat, partitions, trait_names)
    end
end

function DataModel(taxa::AbstractArray{<:AbstractString},
                   data::AbstractMatrix{Float64})
    return DataModel(taxa, data, DEFAULT_TRAIT_NAME)
end

function trait_dimensions(dm::DataModel)
    return [length(x) for x in dm.partitions]
end

