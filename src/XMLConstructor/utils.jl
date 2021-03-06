function df_to_matrix(df::DataFrame) #Converts a data frame with missing values to a matrix with NaN
    taxa = Vector{String}(df[!, :taxon])
    n, p = size(df)
    p = p - 1
    nms = names(df)

    data = fill(NaN, n, p)
    for i = 1:p
        for j = 1:n
            x = df[j, nms[i + 1]]
            if !ismissing(x)
                data[j, i] = x
            end
        end
    end
    return taxa, data
end



function use_dates!(bx::BEASTXMLElement)
    data = get_data(bx)
    newick = get_newick(bx)
    data.use_times = true
    newick.fix_tree = false
end

function add_trait!(bx::BEASTXMLElement, data::Matrix{Float64}, trait::String)
    p = size(data, 2)
    data_el = get_data(bx)
    tm_el = get_treeModel(bx)

    add_trait!(data_el, data, trait)
    add_leaf_param!(tm_el, trait, p)
end

function get_id(el::XMLElement)
    return attribute(el, bn.ID)
end

function set_id!(el::XMLElement, id::String)
    set_attribute(el, bn.ID, id)
end



function add_matrix_parameter(pel::XMLElement, M::AbstractArray{T, 2};
        id::String = "") where T <: Number
    mat_el = new_element(bn.MATRIX_PARAMETER)
    if id != ""
        set_attribute(mat_el, bn.ID, id)
    end
    n, p = size(M)
    for i = 1:n
        param_el = new_element(bn.PARAMETER)
        set_attribute(param_el, bn.VALUE, join(M[i, :], ' '))
        add_child(mat_el, param_el)
    end
    add_child(pel, mat_el)
    return mat_el
end

function add_diagonal_matrix(pel::XMLElement, M::Vector{T};
        id::String = "", lower::String = "") where T <: Number

    mat_el = new_child(pel, bn.DIAGONAL_MATRIX)
    p_el = add_parameter(mat_el, value = M, id=id)
    if lower != ""
        set_attribute(p_el, bn.LOWER, lower)
    end
    return mat_el
end

# function add_parameter(pel::XMLElement, val::Vector{<:Number};
#                        lower::Vector{<:Number} = Float64[])
#     el = new_element(bn.PARAMETER)
#     set_attribute(el, bn.VALUE, join(val, ' '))
#     if length(lower) > 0
#         if length(lower) == 1 || length(lower) == length(val)
#             set_attribute(el, bn.LOWER, join(lower, ' '))
#         else
#             error("incompatible dimension between the parameter values " *
#                 "($(length(val))) and lower bounds ($(length(lower)))")
#         end
#     end

#     add_child(pel, el)
# end


function add_ref_el(pel::XMLElement, el::XMLElement;
            new_name::String = name(el))
    ref_el = reference_element(el, new_name)
    add_child(pel, ref_el)
    return ref_el
end

function add_ref_el(pel::XMLElement, name::String, id::String)
    ref_el = new_child(pel, name)
    set_attribute(ref_el, bn.IDREF, id)
end

function add_ref_els(pel::XMLElement, els::Array{XMLElement})
    for el in els
        add_ref_el(pel, el)
    end
end

function reference_element(el::XMLElement)
    nm = name(el)
    return reference_element(el, nm)
end

function reference_element(el::XMLElement, nm::String)
    id = attribute(el, bn.ID)
    if isnothing(id)
        id = attribute(el, bn.IDREF)
    end

    if isnothing(id)
        @warn "Element does not have id. Adding whole element."
        return el
    else
        ref_el = new_element(nm)
        set_attribute(ref_el, bn.IDREF, id)
        return ref_el
    end
end

function xml_vec(n::Int)
    v = Vector{XMLOrNothing}(undef, n)
    fill!(v, nothing)
    return v
end

function check_posdef(x::AbstractArray{Float64, 2})
    if !isposdef(x)
        error("Matrix must be positive definite.")
    end
end