export LocalOperator, X, Y, Z
using MacroTools

"""
    LocalOperator{T, N, O, MT}

Local operator on `N` qubits/sites with `O` eigen states.
"""
struct LocalOperator{T <: Number, N, O, MT <: AbstractMatrix{T}} <: AbstractMatrix{T}
    val::MT
end

LocalOperator(N::Int, O::Int, val::MT) where {T, MT <: AbstractMatrix{T}} = LocalOperator{T, N, O, MT}(val)

const X = LocalOperator(1, 2, Pauli.X{ComplexF64}())
const Y = LocalOperator(1, 2, Pauli.Y{ComplexF64}())
const Z = LocalOperator(1, 2, Pauli.Z{ComplexF64}())

Base.summary(io::IO, op::LocalOperator{T, N, O}) where {T, N, O} = print(io, "LocalOperator{$T, $N, $O}")

# inherit other array interface
MacroTools.@forward LocalOperator.val Base.getindex,
                                      Base.size, Base.length
                                      Base.stride, Base.strides


Base.eltype(::LocalOperator{T}) where T = T
Base.:(==)(lhs::LocalOperator{T, N, O}, rhs::LocalOperator{T, N, O}) where {T, N, O} = lhs.val == rhs.val

Base.getindex(A::LocalOperator{<:Number, N}, lhs::NTuple{N, <:Integer}, rhs::NTuple{N, <:Integer}) where N =
    getindex(A.val, to_indices(A, (lhs, rhs))...)

@generated function Base.to_index(A::LocalOperator{<:Number, N, O}, i::NTuple{N, <:Integer}) where {N, O}
    ex = :(i[1] - 1)
    for k = 2:N
        factor = O^(k-1)
        ex = :($ex + $factor * (i[$k] - 1))
    end
    :($ex + 1)
end

Base.:(*)(lhs::LocalOperator{T, N, O, MT1}, rhs::LocalOperator{T, N, O, MT2}) where {T, N, O, MT1, MT2} =
    LocalOperator(N, O, lhs.val * rhs.val)
