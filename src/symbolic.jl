export @sum

function sum_expr(ex)
    @match ex begin
        alpha_Number * expr_ => :($alpha * $(sum_expr(expr)))
        expr_ * alpha_Number => :($alpha * $(sum_expr(expr)))

        expr1_ + expr2_ => :($(sum_expr(expr1)) + $(sum_expr(expr2)))
        expr1_ - expr2_ => :($(sum_expr(expr1)) - $(sum_expr(expr2)))
        _ => :(sum($(make_subham(SubHamExpr(ex)))))
    end
end

"""
    @sum local_hamiltonian

Sum over a region defined in `local_hamiltonian`.
"""
macro sum(expr)
    sum_expr(expr)
end

abstract type Operator end

Base.eltype(x::Operator) = eltype(x.expr)

"""
    Sum{H} <: Operator

Symbolic summation.

```math
\\sum H
```
"""
struct Sum{H} <: Operator
    expr::H
end

# NOTE: prod only operates on surround/face?
"""
    Prod{H} <: Operator

Symbolic product.

```math
\\prod H
```
"""
struct Prod{H} <: Operator
    expr::H
end

struct ScalarProd{T <: Number, H} <: Operator
    val::T
    expr::H
end

ScalarProd(val::Number, ex::ScalarProd) = ScalarProd(val * ex.val, ex.expr)

abstract type AddOrSub{LHS, RHS} <: Operator end

struct Add{LHS, RHS} <: AddOrSub{LHS, RHS}
    lhs::LHS
    rhs::RHS
end

struct Sub{LHS, RHS} <: AddOrSub{LHS, RHS}
    lhs::LHS
    rhs::RHS
end

Base.sum(x::SubHam) = Sum(x)
# move scalar out of sum
Base.sum(x::ScalarProd{<:Number}) = ScalarProd(x.val, Sum(x.expr))
# move add/sub out of sum
Base.sum(x::Add) = Add(Sum(x.lhs), Sum(x.rhs))
Base.sum(x::Sub) = Sub(Sum(x.lhs), Sum(x.rhs))

import Base: *, +, -
*(lhs::Number, rhs::Sum) = ScalarProd(lhs, rhs)
*(lhs::Number, rhs::SubHam) = ScalarProd(lhs, rhs)
*(lhs::Number, rhs::ScalarProd) = ScalarProd(lhs * rhs.val, rhs.expr)
*(lhs::Number, rhs::Add) = Add(lhs * rhs.lhs, lhs * rhs.rhs)
*(lhs::Number, rhs::Sub) = Sub(lhs * rhs.lhs, lhs * rhs.rhs)
# scalar multiplication is exchangable
*(lhs::Sum, rhs::Number) = *(rhs, lhs)
*(lhs::SubHam, rhs::Number) = *(rhs, lhs)
*(lhs::AddOrSub, rhs::Number) = *(rhs, lhs)

# minus
-(x::SubHam) = ScalarProd(-1, x)
-(x::Sum) = ScalarProd(-1, x)
-(x::ScalarProd) = ScalarProd(-x.val, x.expr)
-(x::Add) = Sub(-x.lhs, x.rhs)
-(x::Sub) = Add(-x.lhs, x.rhs)

+(lhs::T1, rhs::T2) where {T1 <: Union{Operator, SubHam}, T2 <: Union{Operator, SubHam}} = Add(lhs, rhs)
-(lhs::T1, rhs::T2) where {T1 <: Union{Operator, SubHam}, T2 <: Union{Operator, SubHam}} = Sub(lhs, rhs)
# merge same op
+(lhs::T, rhs::T) where {T <: Union{Operator, SubHam}} = lhs === rhs ? ScalarProd(2, lhs) : Add(lhs, rhs)
-(lhs::T, rhs::T) where {T <: Union{Operator, SubHam}} = lhs === rhs ? ScalarProd(2, lhs) : Sub(lhs, rhs)

# add/sub merge rules

## ScalarProd
+(lhs::ScalarProd{<:Number, T}, rhs::T) where {T <: Union{Operator, SubHam}} = lhs.expr === rhs ? ScalarProd(lhs.val+1, lhs.expr) : Add(lhs, rhs)
-(lhs::ScalarProd{<:Number, T}, rhs::T) where {T <: Union{Operator, SubHam}} = lhs.expr === rhs ? ScalarProd(lhs.val-1, lhs.expr) : Sub(lhs, rhs)
+(lhs::ScalarProd{<:Number, T}, rhs::ScalarProd{<:Number, T}) where T =
    lhs.expr === rhs.expr ? ScalarProd(lhs.val + rhs.val, lhs.expr) : Add(lhs, rhs)
-(lhs::ScalarProd{<:Number, T}, rhs::ScalarProd{<:Number, T}) where T =
    lhs.expr === rhs.expr ? ScalarProd(lhs.val - rhs.val, lhs.expr) : Sub(lhs, rhs)

+(lhs::ScalarProd{<:Number, T}, rhs::Sum{T}) where T =
    lhs.expr === rhs.expr ? ScalarProd(lhs.val + 1, rhs) : Add(lhs, rhs)
-(lhs::ScalarProd{<:Number, T}, rhs::Sum{T}) where T =
    lhs.expr === rhs.expr ? ScalarProd(lhs.val - 1, rhs) : Sub(lhs, rhs)

+(lhs::ScalarProd{<:Number, Sum{T}}, rhs::Add{Sum{T}, <:Operator}) where T =
    lhs.expr === rhs.lhs.expr ? Add(ScalarProd(lhs.val+1, lhs.expr), rhs.rhs) : Add(lhs, rhs)
+(lhs::ScalarProd{<:Number, Sum{T}}, rhs::Add{<:Operator, Sum{T}}) where T =
    lhs.expr === rhs.rhs.expr ? Add(ScalarProd(lhs.val+1, lhs.expr), rhs.lhs) : Add(lhs, rhs)
+(lhs::ScalarProd{<:Number, Sum{T}}, rhs::Sub{Sum{T}, <:Operator}) where T =
    lhs.expr === rhs.lhs.expr ? Sub(ScalarProd(lhs.val+1, lhs.expr), rhs.rhs) : Add(lhs, rhs)
+(lhs::ScalarProd{<:Number, Sum{T}}, rhs::Sub{<:Operator, Sum{T}}) where T =
    lhs.expr === rhs.rhs.expr ? Add(ScalarProd(lhs.val-1, lhs.expr), rhs.lhs) : Add(lhs, rhs)

-(lhs::ScalarProd{<:Number, Sum{T}}, rhs::Add{Sum{T}, <:Operator}) where T =
    lhs.expr === rhs.lhs.expr ? Sub(ScalarProd(lhs.val-1, lhs.expr), rhs.rhs) : Sub(lhs, rhs)
-(lhs::ScalarProd{<:Number, Sum{T}}, rhs::Add{<:Operator, Sum{T}}) where T =
    lhs.expr === rhs.rhs.expr ? Sub(ScalarProd(lhs.val-1, lhs.expr), rhs.lhs) : Sub(lhs, rhs)
-(lhs::ScalarProd{<:Number, Sum{T}}, rhs::Sub{Sum{T}, <:Operator}) where T =
    lhs.expr === rhs.lhs.expr ? Add(ScalarProd(lhs.val-1, lhs.expr), rhs.rhs) : Sub(lhs, rhs)
-(lhs::ScalarProd{<:Number, Sum{T}}, rhs::Sub{<:Operator, Sum{T}}) where T =
    lhs.expr === rhs.rhs.expr ? Sub(ScalarProd(lhs.val+1, lhs.expr), rhs.lhs) : Sub(lhs, rhs)

+(lhs::T, rhs::ScalarProd{<:Number, T}) where T = +(rhs, lhs)
-(lhs::T, rhs::ScalarProd{<:Number, T}) where T = -(rhs, lhs)
+(lhs::Sum, rhs::ScalarProd) = +(rhs, lhs)
+(lhs::Add, rhs::ScalarProd) = +(rhs, lhs)
+(lhs::Sub, rhs::ScalarProd) = +(rhs, lhs)
-(lhs::Add, rhs::ScalarProd) = +(rhs, lhs)
-(lhs::Sub, rhs::ScalarProd) = +(rhs, lhs)


################################# Geometry #####################################

struct OnRegion{E, R} <: Operator
    expr::E
    region::R

    OnRegion(e::E, r::R) where {E, R} = issum(e) ? new{E, R}(e, r) : error("Expect a sum expression, got $e")
end

Base.eltype(x::OnRegion) = eltype(x.expr)

issum(ex::Add) = issum(ex.lhs) && issum(ex.rhs)
issum(ex::Sub) = issum(ex.lhs) && issum(ex.rhs)
issum(ex::Sum) = true
issum(ex::ScalarProd) = issum(ex.expr)
issum(ex::Prod) = false
issum(ex::SubHam) = false

# NOTE: this is a workaround
(ex::Add)(region) = OnRegion(ex, region)
(ex::Sub)(region) = OnRegion(ex, region)
(ex::Sum)(region) = OnRegion(ex, region)
(ex::ScalarProd)(region) = OnRegion(ex, region)

################################# Printing #####################################
Base.show(io::IO, x::Sum{<:SubHam}) = print(io, "sum(", x.expr, ")")
Base.show(io::IO, x::Prod{<:SubHam}) = print(io, "prod(", x.expr, ")")
Base.show(io::IO, x::ScalarProd) = print(io, x.val, "*", x.expr)
Base.show(io::IO, x::Add) = print(io, x.lhs, " + ", x.rhs)
Base.show(io::IO, x::Sub) = print(io, x.lhs, " - ", x.rhs)

function Base.show(io::IO, x::OnRegion)
    printstyled(io, "Expr:\n"; bold=true)
    println(io, x.expr)
    printstyled(io, "Region:\n"; bold=true)
    print(io, x.region)
end

using LaTeXStrings
import LaTeXStrings: latexstring

function _latexstring(x::Sum{<:SubHam})
    index_str = join(x.expr.placeholders, ", ")
    "\\sum_{$(index_str)} $(_latexstring(x.expr))"
end

function _latexstring(x::Prod{<:SubHam})
    index_str = join(x.expr.placeholders, ", ")
    "\\prod_{$(index_str)} $(_latexstring(x.expr))"
end

_latexstring(x::ScalarProd) = "$(x.val) $(_latexstring(x.expr))"
_latexstring(x::Add) = "$(_latexstring(x.lhs) + _latexstring(x.rhs))"
_latexstring(x::Sub) = "$(_latexstring(x.lhs) - _latexstring(x.rhs))"
latexstring(x::Operator) = latexstring("\$\$", _latexstring(x), "\$\$")

Base.show(io::IO, ::MIME"text/latex", x::Operator) = print(io, latexstring(x))

function Base.show(io::IO, ::MIME"text/latex", x::OnRegion)
    printstyled(io, "Expr:\n"; bold=true)
    println(io, latexstring(x.expr))
    printstyled(io, "Region:\n"; bold=true)
    print(io, x.region)
end
