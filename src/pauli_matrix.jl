"""
This module defines types for Pauli Matrix
"""
module Pauli

abstract type PauliMatrix{T <: Complex} <: AbstractMatrix{T} end

struct X{T} <: PauliMatrix{T} end
struct Y{T} <: PauliMatrix{T} end
struct Z{T} <: PauliMatrix{T} end

Base.Matrix(x::PauliMatrix{T}) where T = Matrix{T}(x)
Base.Matrix{T}(::X) where T = T[0 1;1 0]
Base.Matrix{T}(::Y) where T = T[0 -im;im 0]
Base.Matrix{T}(::Z) where T = T[1 0;0 -1]

using SparseArrays

SparseArrays.SparseMatrixCSC(x::PauliMatrix{T}) where T = SparseMatrixCSC{T}(x)
SparseArrays.SparseMatrixCSC{T}(x::X) where T = sparse(Matrix{T}(x))
SparseArrays.SparseMatrixCSC{T}(x::Y) where T = sparse(Matrix{T}(x))
SparseArrays.SparseMatrixCSC{T}(x::Z) where T = sparse(Matrix{T}(x))

SparseArrays.nonzeros(::X{T}) where T = T[1, 1]
SparseArrays.nonzeros(::Y{T}) where T = T[im, -im]
SparseArrays.nonzeros(::Z{T}) where T = T[1, -1]

function Base.checkbounds(::Type{Bool}, x::PauliMatrix, ind1, ind2)
    1 <= ind1 <= 2 || throw(BoundsError(x, (ind1, ind2)))
    1 <= ind2 <= 2 || throw(BoundsError(x, (ind1, ind2)))
    true
end

function Base.getindex(x::X{T}, inds...) where T
    @boundscheck checkbounds(Bool, x, inds...)

    inds == (1, 2) || inds == (2, 1) ? one(T) : zero(T)
end

function Base.getindex(x::Y{T}, inds...) where T
    @boundscheck checkbounds(Bool, x, inds...)

    inds == (1, 2) || inds == (2, 1) ? T(-im) : zero(T)
end

function Base.getindex(x::Z{T}, inds...) where T
    @boundscheck checkbounds(Bool, x, inds...)

    inds == (1, 1) ? -one(T) :
    inds == (2, 2) ? one(T)  :
    zero(T)
end

Base.size(::PauliMatrix) = (2, 2)

Base.summary(io::IO, ::X{T}) where T = print(io, "Pauli X{$T}")
Base.summary(io::IO, ::Y{T}) where T = print(io, "Pauli Y{$T}")
Base.summary(io::IO, ::Z{T}) where T = print(io, "Pauli Z{$T}")

end
