mutable struct ScaledOrthogonalMatrix <: MyXMLElement
    el::XMLOrNothing
    scale::Parameter
    U::MatrixParameter
    id::String

    function ScaledOrthogonalMatrix(scale::Parameter, U::MatrixParameter,
                                    id::String;
                                    tol::Float64=1e-10)

        scale_data = scale.val
        U_data = U.mat
        k = length(scale_data)
        k2, p = size(U_data)

        @assert k == k2

        # check U is orthonormal
        for i = 1:k
            ui = @view U_data[i, :]

            if abs(dot(ui, ui) - 1.0) > tol
                error("matrix is not orthonormal")
            end

            for j = (i + 1):k
                uj = @view U_data[j, :]
                if abs(dot(ui, uj)) > tol
                    error("matrix is not orthonormal")
                end
            end
        end

        return new(nothing, scale, U, id)
    end

end

function ScaledOrthogonalMatrix(L::Matrix{Float64}, id, scale_id, U_id)
    Lsvd = svd(L)
    return ScaledOrthogonalMatrix(Parameter(Lsvd.S, scale_id),
                                  MatrixParameter(Lsvd.Vt, U_id),
                                  id)
end

function make_xml(som::ScaledOrthogonalMatrix)
    make_xml(som.scale)
    make_xml(som.U)
    el = new_element(bn.SCALED_MATRIX)
    set_id!(el, som.id)
    scale_el = new_child(el, bn.SCALE)
    add_ref_el(scale_el, som.scale.el)
    mat_el = new_child(el, bn.MATRIX)
    add_ref_el(mat_el, som.U.el)

    som.el = el
    return el
end

function get_loggables(som::ScaledOrthogonalMatrix)
    make_xml(som)
    return [get_loggables(som.scale); som.el]
end

function set_mat!(som::ScaledOrthogonalMatrix, mat::Matrix{Float64})
    s = svd(mat)
    set_mat!(som.U, s.Vt)
    set_value(som.scale, s.S)
end

function set_scale!(som::ScaledOrthogonalMatrix, scale::Vector{Float64})
    set_value(som.scale, scale)
end

function size(som::ScaledOrthogonalMatrix, args...)
    return size(som.U, args...)
end

function name(::ScaledOrthogonalMatrix)
    return bn.SCALED_MATRIX
end



