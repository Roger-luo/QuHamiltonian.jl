export Sites, sitetype, data
export Bit, Spin, Half
export up, down

################################ Site Label ####################################
"""
    SiteLabel
Site Labels, type of each element, like Spin, Bit, etc.
"""
abstract type SiteLabel end

"""
    BitSiteLabel <: SiteLabel
Binary Site, like Spin, Bit, etc.
"""
abstract type BitSiteLabel <: SiteLabel end

"""
    nstate(label)

Returns number of different states of this label.
"""
nstate(::Type{<:BitSiteLabel}) = 2

#################################################################
# Interfaces

"""
    up(BitSiteLabel) -> tag
up tag for this label. e.g. `1` for `Bit`, `0.5` for `Half`.
"""
function up end

"""
    down(BitSiteLabel) -> tag
down tag for this label. e.g. `0` for `Bit`, `-0.5` for `Half`.
"""
function down end

"""
    Spin <: BitSiteLabel
Binary State, Spin. Has two state: -1, 1
"""
abstract type Spin <: BitSiteLabel end

"""
    Bit <: BitSiteLabel
Binary State, Bit. Has two state: 0, 1
"""
abstract type Bit <: BitSiteLabel end

"""
    Half <: BitSiteLabel
Binary State, Half. Has two state: -0.5, 0.5
"""
abstract type Half <: BitSiteLabel end

Base.eltype(::Type{Spin}) = Float64
up(::Type{Spin}) = convert(eltype(Spin), 1.0f0)
down(::Type{Spin}) = convert(eltype(Spin), -1.0f0)

Base.eltype(::Type{Bit}) = Float64
up(::Type{Bit}) = convert(eltype(Bit), 1.0f0)
down(::Type{Bit}) = convert(eltype(Bit), 0.0f0)

Base.eltype(::Type{Half}) = Float64
up(::Type{Half}) = convert(eltype(Half), 0.5f0)
down(::Type{Half}) = convert(eltype(Half), -0.5f0)

################################## Sites #######################################

"""
    AbstractSites{Label, T, N} <: AbstractArray{T, N}
Sites are arrays with certain labels
"""
abstract type AbstractSites{Label, T, N} <: AbstractArray{T, N} end

"""
    sitetype(sites)

get sitetype from subtypes of `AbstractSites`
"""
function sitetype end

"""
    data(sites)
get data of this `sites`
"""
function data end


# TODO: type promote to normal arrays when mixed with normal arrays.

"""
    Sites{L <: SiteLabel, T, N} <: AbstractSites{L, T, N}

Lattice Sites are Julia Arrays with certain Label, it is able to use array
interface.

# Example

    Sites(Bit, 2, 2)
    Sites(Bit, (2, 2))
    rand(Bit, 2, 2)
"""
mutable struct Sites{L <: SiteLabel, T, N} <: AbstractSites{L, T, N}
    data::Array{T, N}
end

Sites{L}(x::Array{T, N}) where {L, T, N} = Sites{L, T, N}(x)

# Enable type inference
Sites(::Type{L}, data::Array{T, N}) where {L, T, N} = Sites{L, T, N}(data)

Sites(::Type{L}, shape::Int...) where L <: SiteLabel = Sites(L, shape)
Sites(::Type{L}, shape::Tuple) where L <: SiteLabel =
    Sites{L, eltype(L), length(shape)}(fill(down(L), shape))

data(s::Sites) = s.data

# use array interface
Base.eltype(x::Sites{L, T}) where {L, T} = T
Base.length(x::Sites) = length(x.data)
Base.ndims(x::Sites) = ndims(x.data)
Base.size(x::Sites) = size(x.data)
Base.size(x::Sites, n::Integer) = size(x.data, n)
Base.axes(x::Sites) = axes(x.data)
Base.axes(x::Sites, d::Integer) = axes(x.data, d)
Base.eachindex(x::Sites) = eachindex(x.data)
Base.stride(x::Sites, k::Integer) = stride(x.data, k)
Base.strides(x::Sites) = strides(x.data)
Base.getindex(x::Sites, index::Integer...) = getindex(x.data, index...)
Base.getindex(x::Sites, index::NTuple{N, T}) where {N, T <: Integer} = getindex(x.data, index...)
Base.setindex!(x::Sites, val, index::Integer...) = setindex!(x.data, val, index...)
Base.setindex!(x::Sites, val, index::NTuple{N, T}) where {N, T <: Integer} = setindex!(x.data, val, index...)

# TODO: value check
Base.getindex(::Type{Bit}, xs...) = Sites{Bit}([xs...])
Base.getindex(::Type{Spin}, xs...) = Sites{Spin}([xs...])
Base.getindex(::Type{Half}, xs...) = Sites{Spin}([xs...])

@inline function Base.copy(b::Sites{L, T, N}) where {L, T, N}
    return Sites{L, T, N}(copy(b.data))
end

sitetype(::Type{Sites{L}}) where L = L
sitetype(x::Sites{L}) where L = L

export reset!, set!

@inline function reset!(b::Sites{L, T, N}) where {L, T, N}
    fill!(b.data, down(L))
    return b
end

@inline function set!(s::Sites{Bit, T, N}, n::Int) where {T, N}
    for i in eachindex(s)
        @inbounds s[i] = ((n>>(i-1)) & 0x1)
    end
    return s
end

@inline function set!(s::Sites{Spin, T, N}, n::Int) where {T, N}
    for i in eachindex(s)
        @inbounds s[i] = 2 * ((n>>(i-1)) & 0x1) - 1
    end
    return s
end

@inline function set!(s::Sites{Half, T, N}, n::Int) where {T, N}
    for i in eachindex(s)
        @inbounds s[i] = ((n>>(i-1)) & 0x1) - 0.5
    end
    return s
end

#######################
# randomize interface
#######################

using Random

function Random.rand!(rng::AbstractRNG, b::Sites{L, T, N}) where {L, T, N}
    rand!(rng, b.data, [up(L), down(L)])
    return b
end

Random.rand!(b::Sites{L, T, N}) where {L <: SiteLabel, T, N} = rand!(GLOBAL_RNG, b)
Random.rand(rng::AbstractRNG, ::Type{L}, shape::Dims) where {L <: SiteLabel} =
    rand!(rng, Sites(L, shape))

####################

@inline function flip!(b::Sites{L, T, N}, index::Integer...) where {L, T, N}
    if b[index...] == up(L)
        b[index...] = down(L)
    else
        b[index...] = up(L)
    end
    return b
end

@inline function randflip!(b::Sites{L, T, N}) where {L, T, N}
    offset = rand(1:length(b))
    return flip!(b, offset)
end

# carry bit
function carrybit!(a::Sites{L}) where L
    @inbounds for i in eachindex(a)
        if a[i] == up(L)
            a[i] = down(L)
        else
            a[i] = up(L)
            break
        end
    end
    a
end

"""
    ndowns(sites)

count the number of down tags, e.g number of `0`s
"""
function ndowns(s::Sites{L}) where {L <:BitSiteLabel}
    count = 0
    @inbounds for i in eachindex(s)
        if s[i] == down(L)
            count += 1
        end
    end
    count
end

function nups(s::Sites{L}) where {L <:BitSiteLabel}
    count = 0
    @inbounds for i in eachindex(s)
        if s[i] == up(L)
            count += 1
        end
    end
    count
end

function n_downs_ups(s::Sites{L}) where {L <:BitSiteLabel}
    count_down = 0
    count_up   = 0
    @inbounds for i in eachindex(s)
        s[i] == up(L) ? (count_up += 1) : (count_down += 1)
    end
    count_down, count_up
end

# TODO: whether this should have side effect?
# will side effect affects performance (extra instance)?
function Base.:(<<)(a::Sites, b::Int)
    for i = 1:b
        carrybit!(a)
    end
    return a
end

##############
# Conversions
##############

import Base: convert

function convert(::Type{Integer}, x::Sites{L, T, N}) where {L, T, N} end

for IntType in (:Int8, Int16, Int32, Int64, Int128, BigInt)
@eval begin
        function convert(::Type{$IntType}, x::Sites{Bit, T, N}) where {T, N}
            if sizeof($IntType) * 8 < length(x)
                throw(Compat.InexactError(:convert, $IntType, x))
            end

            sum($IntType(each) << (i-1) for (i, each) in enumerate(x))
        end

        function convert(::Type{$IntType}, x::Sites{Spin, T, N}) where {T, N}
            if sizeof($IntType) * 8 < length(x)
                throw(Compat.InexactError(:convert, $IntType, x))
            end

            sum($IntType(div(each+1, 2)) << (i-1) for (i, each) in enumerate(x))
        end

        function convert(::Type{$IntType}, x::Sites{Half, T, N}) where {T, N}
            if sizeof($IntType) * 8 < length(x)
                throw(Compat.InexactError(:convert, $IntType, x))
            end

            sum($IntType(each+0.5) << (i-1) for (i, each) in enumerate(x))
        end
    end
end
