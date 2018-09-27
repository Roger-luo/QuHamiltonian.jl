export @subham

using MacroTools
using MacroTools: postwalk, @match

# Expr(:placeholder, :i, :j)
function is_hamiltonian_expr(expr)
    @match expr begin
        # Z[i]
        # A[i, j, k]
        _[__] => true

        # Z[i]Z[j]
        # A[i, j, k]B[i, j]
        *(args__) => all(map(is_hamiltonian_expr, args))
        ⊗(args__) => all(map(is_hamiltonian_expr, args))

        # Z[i]Z[j] where i where j
        ex_ where placeholder_ => is_hamiltonian_expr(ex)
        ex_ where placeholders__ => is_hamiltonian_expr(ex)

        # Expr(:call, [])
        _ => false
    end
end

is_hamiltonian_expr(:(Z[i]Z[j] where {i, j <: nearest}))

@capture(:(Z[i]Z[j] where {i, j <: nearest}), ex_ where placeholder__)

placeholder

function ref_list(expr)
    is_hamiltonian_expr(expr) || throw(Meta.ParseError("Invalid Syntax, expect hamiltonian expr, got $expr"))
    ref_list = Expr[]
    postwalk(expr) do x

    end
end

struct PlaceHolderExpr{N}
    head::Symbol
    names::NTuple{N, Symbol}
end

PlaceHolderExpr(name::Vector{Symbol}) = PlaceHolderExpr(:free, name)
PlaceHolderExpr(name::Symbol) = PlaceHolderExpr(:free, [name])

struct SubHamExpr{N}
    placeholders::PlaceHolderExpr{N}
    val::Vector{Any}
end

function SubHamExpr(ex::Expr)
    @match ex begin
        # Z[i, j]
        name_[vars__] => SubHamExpr(PlaceHolderExpr(vars))

        *(args__) => SubHamExpr(args)
    end
end

function SubHamExpr(ex::Vector)
    args = []
    for each in ex
        push!(args, SubHamExpr(ex))
    end
    SubHamExpr()
end
