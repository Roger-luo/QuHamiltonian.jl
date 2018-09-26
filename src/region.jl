export UndefRegion
import Lattices

"""
    UndefRegion{N}

Undefined region on `N` sites.
"""
struct UndefRegion{N} end
UndefRegion(n::Int) = UndefRegion{n}()

struct UndefRegionIter{N}
    UndefRegionIter(x::UndefRegion{N}) where N = new{N}()
end

Lattices.sites(x::UndefRegion) = UndefRegionIter(x)
Lattices.edges(::UndefRegion; length=1) = error("There is no edges for UndefRegion Region")
Base.iterate(it::UndefRegionIter{N}) where N = N == 0 ? nothing : (1, 1)
Base.iterate(it::UndefRegionIter{N}, i) where N = i == N + 1 ? nothing : (i + 1, i + 1)

Base.getindex(it::UndefRegionIter{N}, i::Tuple{Int}) where N = getindex(it, first(i))
Base.getindex(it::UndefRegionIter{N}, i::Int) where N = mod1(i, N)
Base.length(it::UndefRegionIter{N}) where N = N
Base.eltype(it::UndefRegionIter) = Int

@generated function Base.collect(::Type{T}, src::UndefRegionIter{N}) where {T, N}
    :(T[i for i in 1:$N])
end

Base.collect(src::UndefRegionIter) = collect(Int, src)
