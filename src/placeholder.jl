export @region, @neighbor, @nearest, @surround, @face
using MacroTools
using MacroTools: @match, postwalk

abstract type PlaceHolder end
name(ph::PlaceHolder) = ph.name
Base.show(io::IO, ph::PlaceHolder) = print(io, name(ph))

const PlaceHolderMap = Dict{Symbol, Type}()

function _placeholder_expr(typename)
    quote
        Core.@__doc__ struct $(typename) <: PlaceHolder
            name::Symbol
        end
    end
end

function _placeholder_var_expr(type_::Type, var)
    :(Core.@__doc__ $(esc(var)) = $(type_(var)))
end

declare_placeholder(type_::Type, var::Symbol) = _placeholder_var_expr(type_, var)

function declare_placeholder(type_::Type, vars)
    @assert vars.head == :tuple || vars.head == :escape "invalid syntax, expect tuple, got $vars ($vars.head)"

    ex = quote end
    for each in vars.args
        def = _placeholder_var_expr(type_, each)
        push!(ex.args, def)
    end
    push!(ex.args, esc(vars))
    ex
end

make_region(alias::Symbol, vars) = declare_placeholder(PlaceHolderMap[alias], vars)
make_neighbor(k::Int, vars) = declare_placeholder(Neighbor{k}, vars)

function register_holder_macro_entry(name::QuoteNode, val)
    mname = Symbol("@$(name.value)")
    docstr =
"""
    @$(name.value) vars...

Declare place holders for `$val`. See doc of `$val` for details.
"""

    quote
        PlaceHolderMap[$name] = $(esc(val))

        export $mname # export this macro
        macro $(esc(name.value))(vars)
            make_region($name, vars)
        end

        @doc $docstr $mname
    end
end

"""
    @register_holder alias => typename
    @register_holder typename alias_list

Define place holder types and register their alias in `PlaceHolderMap`.

### Example

register a simple `Vertex` type with alias `vertex`.

```julia
@register_holder :vertex => Vertex
```

register a parametric type `Neighbor{K}` with different alias.

```julia
@register_holder Neighbor{K} begin
    :nearest => Neighbor{1}
    :next_nearest => Neighbor{2}
end
```
"""
macro register_holder end

macro register_holder(typename, list)
    ex = _placeholder_expr(typename)

    for each in list.args
        if @capture(each, name_ => val_)
            def = register_holder_macro_entry(name, val)
            push!(ex.args, def)
        end
    end
    # push expr return value
    push!(ex.args, true)
    ex
end

macro register_holder(pair)
    @capture(pair, name_ => typename_)
    quote
        $(_placeholder_expr(typename))

        $(register_holder_macro_entry(name, typename))
        true
    end
end

################################################################################

"""
    Vertex <: PlaceHolder

Place holder for vertex on lattice. This represents
all the vertexes on a lattice/region.
"""
@register_holder :vertex => Vertex

"""
    Neighbor{K} <: PlaceHolder

Place holder for near-neighbor correlation. This represents
k-th near-neighbor on a lattice/region. e.g `Neighbor{1}` means
the nearest neighbor.
"""
@register_holder Neighbor{K} begin
    :nearest => Neighbor{1}
    :next_nearest => Neighbor{2}
end

"""
    Surround <: PlaceHolder

Place holder for surrounding area of a vertex.
"""
@register_holder :surround => Surround
@register_holder :face => Face

################################################################################
# Factory Macros
macro region(alias::QuoteNode, vars)
    make_region(alias.value, vars)
end

macro neighbor(k::Int, vars)
    make_neighbor(k, vars)
end
