struct NumberOfBody{K} end
NumberOfBody(::NumberOfBody{K}) where K = NumberOfBody{K}()
NumberOfBody(::LocalOperator{T, N}) where {T, N} = NumberOfBody{N}()
NumberOfBody(x::SubHam) = NumberOfBody{length(x.placeholders)}()
NumberOfBody(x::Operator) = NumberOfBody(x.expr)

struct PlaceHolderType{T} end
PlaceHolderType(::SubHam{PT}) where PT = PlaceHolderType{PT}()

# NOTE: only operator acts on a single SubHam
PlaceHolderType(x::Operator) = PlaceHolderType(x.expr)
# PlaceHolderType(x::AddOrSub) =

"""
    LOpEigen{K}

Number of eigen states in local operator.
"""
struct LOpEigen{K} end

LOpEigen(::LocalOperator{T, N, O}) where {T, N, O} = LOpEigen{O}()
LOpEigen(x::Operator) = LOpEigen(x.expr)
LOpEigen(x::AddOrSub) = (LOpEigen(x.lhs); LOpEigen(x.rhs))
LOpEigen(x::SubHam{<:Tuple, <:LocalOperator}) = LOpEigen(x.val)

function LOpEigen(x::SubHam{<:Tuple, <:Tuple})
    n = LOpEigen(first(x.val))

    # check if all the local operator has the same number of eigens
    for each in x.val
        n === LOpEigen(each) || error("number of eigen states in given local operators are not the same")
    end
    n
end

struct LocalOperatorValType{T, OpType} end
LocalOperatorValType(::MT) where {T, MT <: AbstractMatrix{T}} = LocalOperatorValType{T, MT}()
LocalOperatorValType(x::LocalOperator) = LocalOperatorValType(x.val)
LocalOperatorValType(x::SubHam{<:Tuple, <:LocalOperator}) = LocalOperatorValType(x.val)
LocalOperatorValType(x::Operator) = LocalOperatorValType(x.expr)

# multiple operators
LocalOperatorValType(x::SubHam{<:Tuple, <:Tuple}) = LocalOperatorValType(x.val...)
LocalOperatorValType(a::SubHam, x::SubHam...) = LocalOperatorValType(LocalOperatorValType(a), LocalOperatorValType(x...))

LocalOperatorValType(x::LocalOperatorValType{T}...) where T = LocalOperatorValType{T, Tuple{_op_type(x...)...}}()
_op_type(::Tuple{}) = ()
_op_type(::LocalOperatorValType{T, Op}) where {T, Op} = Op
_op_type(::LocalOperatorValType{T, Op1}, ::LocalOperatorValType{T, Op2}, ops...) where {T, Op1, Op2} = (Op1, Op2, _op_type(ops)...)
