mat(ex::SubHam{<:Tuple, <:LocalOperator}) = ex.val

# convert config to indices
Base.to_indices(::Type{L}, inds...) where L <: SiteLabel = map(x->to_indices(L, x), inds)

Base.to_indices(::Type{Bit}, val::Number) = Int(val + 1)
Base.to_indices(::Type{Spin}, val::Number) = Int(div(val + 1, 2))
Base.to_indices(::Type{Half}, val::Number) = Int(val + 0.5)

@generated function Base.to_indices(::Type{Bit}, val::NTuple{N, T}) where {N, T}
    ex = :(val[1])
    for i = 2:N
        factor = 2^(i-1)
        ex = :(val[$i] * $factor + $ex)
    end
    :(ex)
end

@generated function Base.to_indices(::Type{Spin}, val::NTuple{N, T}) where {N, T}
    ex = :(div(val[1] + 1, 2))
    for i = 2:N
        factor = 2^(i-1)
        ex = :(div(val[$i] + 1, 2) * $factor + $ex)
    end
    :(Int($ex))
end

@generated function Base.to_indices(::Type{Half}, val::NTuple{N, T}) where {N, T}
    ex = :(val[1] + 0.5)
    for i = 2:L
        factor = 2^(i-1)
        ex = :((val[$i] + 0.5) * $factor + $ex)
    end
    :(Int($ex))
end

function is_other_sites_the_same(::NumberOfBody{1}, site, region, lhs::Sites, rhs::Sites)
    @inbounds for other_site in Lattices.sites(region)
        if other_site != site
            s = region[other_site]
            lhs[s] == rhs[s] || return false
        end
    end
    true
end

function is_other_sites_the_same(::NumberOfBody, sites, region, lhs::Sites, rhs::Sites)
    for other_site in Lattices.sites(region)
        # NOTE: we seperate single site because of the following:
        #
        # # chain lattice
        # julia> 1 in 1
        #        true
        #
        # # square lattice
        # julia> (1, 1) in (1, 1)
        #        false
        if !(other_site in sites)
            s = region[other_site]
            lhs[s] == rhs[s] || return false
        end
    end
    true
end


Base.getindex(ex::OnRegion, lhs::Sites, rhs::Sites) = eval_expr(ex, lhs, rhs)

## eval Operator expressions
# unzip OnRegion
eval_expr(ex::OnRegion, lhs::Sites, rhs::Sites) =
    eval_expr(ex.expr, ex.region, lhs, rhs)

# 1. dispatch to specific local operator
eval_expr(ex::Operator, region, lhs::Sites, rhs::Sites) =
    eval_expr(ex, region, LocalOperatorValType(ex), PlaceHolderType(ex), lhs, rhs)

# 2. if there is no specific method for the operator, then fall back to general methods
eval_expr(ex::Operator, region, ::LocalOperatorValType, ::PlaceHolderType, lhs, rhs) =
    eval_expr(ex, region, LOpEigen(ex), NumberOfBody(ex), PlaceHolderType(ex), lhs, rhs)

# H1 + H2
eval_expr(ex::Add, region, lhs::Sites, rhs::Sites) =
    eval_expr(ex.lhs, region, lhs, rhs) + eval_expr(ex.rhs, region, lhs, rhs)
# H1 - H2
eval_expr(ex::Sub, region, lhs::Sites, rhs::Sites) =
    eval_expr(ex.lhs, region, lhs, rhs) - eval_expr(ex.rhs, region, lhs, rhs)

# α * H
eval_expr(ex::ScalarProd, region, lhs::Sites, rhs::Sites) =
    ex.val * eval_expr(ex.expr, region, lhs, rhs)

# sum(A[i]) on vertex
# NOTE:
# The summation looks like
#
# U x I x I + I x U x I + I x I x U
#
# which means the 2x2 local operator U
# defines which config should be the same/or not
#
# So we will check if position (excludes U) are the same in each loop
# this determine if this block is zero.
# TODO:
# 1. for diagnol local operators, we can directly check if all configurations
# are the same.
# 2. for local operators like X, we can directly check if all configurations
# are opposite (regarding to LOpEigen{2} operators)
function eval_expr(
    expr::Sum{<:SubHam}, region,
    ::LOpEigen{2}, n::NumberOfBody{1}, # 2x2 local operator
    ::PlaceHolderType{Tuple{Vertex}},
    lhs::Sites{L}, rhs::Sites{L}) where L

    Aij = zero(eltype(expr))

    @inbounds for site in Lattices.sites(region)
        s = region[site]
        lhs_ind = to_indices(L, lhs[s])
        rhs_ind = to_indices(L, rhs[s])

        # check identities
        Aij += is_other_sites_the_same(n, site, region, lhs, rhs) ? expr.expr[lhs_ind, rhs_ind] : 0
    end
    Aij
end

######################### specialized methods ##################################
# # sum(X[i]) over all vertex
# function eval_expr(
#     ex::Sum{<:SubHam}, region,
#     ::LocalOperatorValType{T, <:Pauli.X},
#     ::PlaceHolderType{Tuple{Vertex}},
#     lhs::Sites{L}, rhs::Sites{L}) where {T, L}
#
#     # check if all the configurations are opposite
#     @inbounds for site in Lattices.sites(region)
#         s = region[site]
#         # NOTE: since this is for Pauli.X
#         # we can assume there is only two
#         # states, so just use !=
#         lhs[s] != rhs[s] || return zero(T)
#     end
#
#     # T(length(lhs))
#     # TODO:
#     # NOT DONE YET
# end


# sum(Z[i]) over all vertex
function eval_expr(
    ex::Sum{<:SubHam}, region,
    ::LocalOperatorValType{T, <:Pauli.Z},
    ::PlaceHolderType{Tuple{Vertex}},
    lhs::Sites{L}, rhs::Sites{L}) where {T, L}

    # check if all the configurations are the same
    @inbounds for site in Lattices.sites(region)
        s = region[site]
        lhs[s] == rhs[s] || return zero(T)
    end

    # down down -> -1
    # up   up   -> +1
    downs, ups = n_downs_ups(lhs)
    T(ups - downs)
end

# sum(Z[i]Z[j]) on neighbors
function eval_expr(
    expr::Sum{<:SubHam}, region,
    ::LOpEigen{2}, n::NumberOfBody{2}, # 2x2 local operator
    ::PlaceHolderType{Tuple{Neighbor{K}, Neighbor{K}}},
    lhs::Sites{L}, rhs::Sites{L}) where {K, L}

    Aij = zero(eltype(expr))
    @inbounds for (site_a, site_b) in Lattices.edges(region, length=K)
        s_a = region[site_a]
        s_b = region[site_b]

        lhs_ind = to_indices(L, lhs[s_a], lhs[s_b])
        rhs_ind = to_indices(L, rhs[s_a], rhs[s_b])

        Aij += is_other_sites_the_same(n, (site_a, site_b), region, lhs, rhs) ? expr.expr[lhs_ind, rhs_ind] : 0
    end
    Aij
end

# partial indexing
# TODO: this should use partial eval
# TODO: add assertion
Base.getindex(ex::OnRegion, lhs::Sites, colon::Colon) = eval_expr(ex.expr, ex.region, lhs, colon)
Base.getindex(ex::OnRegion, colon::Colon, rhs::Sites) = eval_expr(ex.expr, ex.region, colon, rhs)

function eval_expr(
    ex, region,
    lhs::Sites{L, T, N}, ::Colon) where {L, T, N}

    # TODO:
    # wrap this as an iterator?
    # for rhs in IterName
    # blabla
    # end

    nelem = 1<<length(region)
    Ai = Vector{eltype(ex)}(undef, nelem)
    rhs = Sites{L, T, N}(fill(down(L), size(lhs)))
    for i = 1:nelem
        Ai[i] = eval_expr(ex, region, lhs, rhs)
        rhs << 1
    end

    Ai
end

function eval_expr(
    ex, region,
    ::Colon, rhs::Sites{L, T, N}) where {L, T, N}

    # TODO:
    # wrap this as an iterator?
    # for rhs in IterName
    # blabla
    # end

    nelem = 1<<length(region)
    Ai = Vector{eltype(ex)}(undef, nelem)
    lhs = Sites{L, T, N}(fill(down(L), size(rhs)))
    for i = 1:nelem
        Ai[i] = eval_expr(ex, region, lhs, rhs)
        lhs << 1
    end

    Ai
end

################################################################################

# TODO:
# sum(A1[i]A2[j]...) on vertex
LocalOperator(x::SubHam{<:Tuple, <:LocalOperator}) = x.val
LocalOperator{T, 1}(x::SubHam{<:Tuple, <:LocalOperator{T, 1}}) where T = x.val
LocalOperator{T, 1}(x::SubHam{<:Tuple, <:Tuple}) where T = prod(LocalOperator{T, 1}, x.val)

# TODO: use a repeated product instead
# 1. A1[i]A2[j] does not go through (1, 1), (2, 2), ..., this should use @subham A1[i]A2[i]
# 2. parse A1[i]A2[i] => (A1 * A2)[i] in make_subham, not in SubHamExpr, we keep A1[i]A2[i] for printing

function eval_expr(ex::Sum{<:SubHam}, region,
    ::LOpEigen{2}, n::NumberOfBody{N}, # 2^Nx2^N local operator
    ::PlaceHolderType{NTuple{N, Vertex}},
    lhs::Sites{L}, rhs::Sites{L}) where {N, L}

    error("Summation over multiple vertex is not supported yet. See issue #1")
    Aij = zero(eltype(ex))

    it = Lattices.sites(region)

    it = Iterators.product((it for i=1:N)...)

    local_h = LocalOperator{eltype(ex), 1}(ex.expr)

    @inbounds for _sites in it
        lhs_ind = Tuple(to_indices(L, lhs[region[site]]) for site in _sites)
        rhs_ind = Tuple(to_indices(L, rhs[region[site]]) for site in _sites)

        Aij += is_other_sites_the_same(n, _sites, region, lhs, rhs) ? ex.expr[lhs_ind, rhs_ind] : 0
    end
    Aij
end

# struct NonZeroIterator
# end
#
# function eval_expr(
#     ex::Sum{<:SubHam}, region,
#     ::LOpEigen{2}, ::NumberOfBody{2}, ::PlaceHolderType{Tuple{Vertex}},
#     lhs::Sites{L}) where L
#
#     Ai = []
#     @inbounds for site in Lattices.sites(region)
#         s = region[site]
#         lhs_ind = to_indices(L, lhs[s])
#         rhs_ind = to_indices(L, rhs[s])
#
#         Aij += expr.expr[lhs_ind, rhs_ind]
#     end
# end
