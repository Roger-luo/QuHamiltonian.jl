export @subham
using MacroTools
using MacroTools: postwalk, @match

merge_placeholder_expr(::Tuple{}, holders::NTuple{N, Symbol}) where N = holders

function merge_placeholder_expr(a::NTuple{N1, Symbol}, b::NTuple{N2, Symbol}) where {N1, N2}
    (a..., Tuple(each for each in b if !(each in a))...)
end

function is_hamiltonian_expr(expr)
    @match expr begin
        # Z[i]
        # A[i, j, k]
        _[__] => true

        # Z[i]Z[j]
        # A[i, j, k]B[i, j]
        *(args__) => all(map(is_hamiltonian_expr, args))
        ⊗(args__) => all(map(is_hamiltonian_expr, args))

        # Expr(:call, [])
        _ => false
    end
end

function decompose(expr)
    is_hamiltonian_expr(expr) || throw(Meta.ParseError("Invalid Syntax, expect hamiltonian expr, got $expr"))

    ref_list = Expr[]
    postwalk(expr) do x
        Meta.isexpr(x, :ref) && push!(ref_list, x)
        x
    end

    ref_list
end

abstract type HamiltonianExpr end

struct SubHamExpr{PT <: Tuple, VT <: Union{Tuple, Symbol}} <: HamiltonianExpr
    placeholders::PT
    val::VT
end

function SubHamExpr(ex::Expr) # A[i]
    @capture(ex, name_[vars__]) && return SubHamExpr(Tuple(vars), name)

    vals = Tuple(SubHamExpr(each) for each in decompose(ex))

    placeholders = ()
    for each in vals
        placeholders = merge_placeholder_expr(placeholders, each.placeholders)
    end
    SubHamExpr(placeholders, vals)
end

macro subham_lower(ex)
    SubHamExpr(ex)
end

function Base.show(io::IO, ex::SubHamExpr{<:Tuple, Symbol})
    print(io, ex.val, "[", join(ex.placeholders, ", "), "]")
end

function Base.show(io::IO, ex::SubHamExpr{<:Tuple, <:Tuple})
    for each in ex.val
        show(io, each)
    end
end

struct SubHam{PT <: Tuple, VT}
    placeholders::PT
    val::VT
    name::Symbol

    SubHam(ph::PT, val::VT) where {PT, VT} = new{PT, VT}(ph, val)
    SubHam(ph::PT, val::VT, name::Symbol) where {PT, VT} = new{PT, VT}(ph, val, name)
end

Base.eltype(x::SubHam{<:Tuple, <:LocalOperator}) = eltype(x.val)
Base.eltype(x::SubHam{<:Tuple, <:Tuple}) = promote_type((eltype(each) for each in x.val)...)

function make_subham(ex::SubHamExpr{<:Tuple, Symbol})
    placeholders = Expr(:tuple, (esc(each) for each in ex.placeholders)...)
    quote
        SubHam($(placeholders), $(esc(ex.val)), $(QuoteNode(ex.val)))
    end
end

function make_subham(ex::SubHamExpr{<:Tuple, <:Tuple})
    placeholders = Expr(:tuple, (esc(each) for each in ex.placeholders)...)
    vals = Expr(:tuple, (make_subham(each) for each in ex.val)...)
    quote
        SubHam($(placeholders), $(vals))
    end
end

macro subham(expr)
    ex = SubHamExpr(expr)
    make_subham(ex)
end

function _plain(h::SubHam{<:Tuple, <:LocalOperator})
    index_str = join(h.placeholders, ", ")
    "$(h.name)[$index_str]"
end

function _plain(h::SubHam{<:Tuple, <:Tuple})
    join(_plain(each) for each in h.val)
end

Base.show(io::IO, h::SubHam) = print(io, _plain(h))

import LaTeXStrings: latexstring

function _latexstring(h::SubHam{<:Tuple, <:LocalOperator})
    index_str = join(h.placeholders, ", ")
    "$(h.name)_{$index_str}"
end

function _latexstring(h::SubHam{<:Tuple, <:Tuple})
    join(_latexstring(each) for each in h.val)
end

latexstring(h::SubHam) = latexstring(_latexstring(h))

Base.show(io::IO, ::MIME"text/latex", h::SubHam) = print(io, latexstring(h))

nplaceholders(h::SubHam) = length(h.placeholders)
placeholder_type(h::SubHam{PT}) where PT = PT

function merge_placeholder(a::Tuple, b::Tuple)
    (a..., Tuple(each for each in b if !(each in a))...)
end

import Base: *, ==

# NOTE: we only check if placeholders are the same type
==(lhs::SubHam, rhs::SubHam) = false

function ==(lhs::SubHam{PT, <:Tuple}, rhs::SubHam{PT, <:Tuple}) where {PT <: Tuple}
    all(map(==, lhs.val, rhs.val))
end

function ==(lhs::SubHam{PT, <:LocalOperator}, rhs::SubHam{PT, <:LocalOperator}) where {PT <: Tuple}
    lhs.val == rhs.val
end

function *(lhs::SubHam, rhs::SubHam)
    placeholders = merge_placeholder(lhs.placeholders, rhs.placeholders)
    SubHam(placeholders, (lhs, rhs))
end

# error
function *(lhs::SubHam{NTuple{2, Neighbor{K}}}, rhs::SubHam) where K
    error("Neighbor correlation can only be applied on 2-local Hamiltonian")
end

# indexing
# single site
Base.getindex(h::SubHam{<:Tuple, <:LocalOperator}, inds...) = getindex(h.val, inds...)

struct IndexTable{N, PT <: Tuple}
    syms::PT
    lhs_inds::NTuple{N, Int}
    rhs_inds::NTuple{N, Int}
end

IndexTable(h::SubHam, lhs, rhs) = IndexTable(h.placeholders, lhs, rhs)

function pick(table::IndexTable, inds::NTuple{K}) where K
    indices = findall(x->in(x, inds), table.syms)
    table.lhs_inds[indices], table.rhs_inds[indices]
end

# TODO: unroll the loop with generated & Cassette
# NOTE: this should directly generate the following, since we have almost everything
# in type.
# A[(i1, i2, ..., iN), (j1, j2, ..., jN)] := Z[i1, j1] * Z[i2, j2] * C[i3, i4, j3, j4] ...
function Base.getindex(h::SubHam{<:Tuple, <:Tuple}, lhs::NTuple{N, Int}, rhs::NTuple{N, Int}) where N
    table = IndexTable(h, lhs, rhs)
    elem = one(eltype(h))
    for each in h.val
        each_lhs, each_rhs = pick(table, each.placeholders)
        elem *= getindex(each, each_lhs, each_rhs)
    end
    elem
end
