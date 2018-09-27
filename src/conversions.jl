Base.Matrix(ex::OnRegion) = Matrix{eltype(ex)}(ex)
Base.Matrix{T}(ex::OnRegion) where T = Matrix{T}(ex, LOpEigen(ex))

function Base.Matrix{T}(expr::OnRegion, ::LOpEigen{2}) where T
    lhs = Sites(Bit, size(expr.region))
    rhs = Sites(Bit, size(expr.region))
    N = length(expr.region)
    H = Matrix{T}(undef, 1<<N, 1<<N)

    @inbounds for i = 1:1<<N
        for j = 1:1<<N
            H[j, i] = expr[lhs, rhs]
            lhs << 1
        end
        rhs << 1
    end
    H
end

using SparseArrays
# TODO: use partial eval instead
# TODO: And parallel the loop
SparseArrays.sparse(ex::OnRegion) = sparse(Matrix(ex))
