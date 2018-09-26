# TODO: finish this!
struct RepeatedProductIterator{N, I}
    iterator::I
end

Iterators.product(N::Int, it) = RepeatedProductIterator(N, it)

Base.IteratorSize(::Type{RepeatedProductIterator{N, I}}) where {N, I} = Base.IteratorSize(I)
