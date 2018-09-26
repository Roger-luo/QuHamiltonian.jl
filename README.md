# Hamiltonian Expression Parser

I was drunk on my birthday, then I wrote this.

This package implements a Domain Specific Language (DSL) for Hamiltonian construction.

## A glance of the interface

```julia
# define some scalar parameters
J = 0.1
h = 0.1

# define place holders
@nearest i, j
@vertex l

H = J * @sum(Z[i]Z[j]) - h * @sum(X[l])
```

Then the `H` object has already contain the whole information of a **TFIM Hamiltonian**
but to use it you need to define a region that provides how placeholders
`i, j` and `l` is evaluated. This can be a geometry information like lattice,
or any other region type that defines the [Region Interface]().

We provide an undefined region type `UndefRegion{N}` for models like **SYK Model**.


Here we just simply use a `Chain` lattice as example

```julia
lattice = Chain(4)
TFI_on_chain = H(lattice) # here we go!
```

the `TFI_on_chain` can be used in various ways since it is actually just a type
contains the whole information about a **TFIM** on chain lattice. e.g

### Access its matrix form

We provide `Sites` type to define a site configuration on lattices, it is actually
just a wrapper of builtin `Array`, but it add an tag about the configuration, e.g
`Bit` for **0, 1** configuration, `Spin` for **-1, 1** configuration and `Half` for
`-0.5, 0.5` configuration.

This will make conversions and processing of configurations much easier. Therefore
you will be able to convert an instance of `Sites` to a `Int`.


```julia
lhs = Bit[0, 0, 1, 0]
rhs = Bit[1, 0, 1, 1]

TFI_on_chain[lhs, rhs] # This gives you the element of given configuration.
```

And of course, convert it to a matrix!

```julia
Matrix(TFI_on_chain)
SparseMatrixCSC(TFI_on_chain)
```

Even GPU! (You need to buy one first.)

```julia
cu(TFI_on_chain) # this will just add a tag of device
CuArray(TFI_on_chain)
```

## Is that all?! Nope.

### Symbolic Computing

Simple symbolic computing is actually supported. So
if you defined some weird Hamiltonian, you can leave
it to the parser and have a beer!

e.g

```julia-repl
julia> @sum(Z[@vertex i]) + 2 * @sum(Z[@vertex i])
3*sum(Z[i])
```

### Interact with other libraries

No hurries, it's coming!

- [ ] conversion to MPO (waiting for iTensor?!)


julia> h = @subham Z[@nearest i]Z[@nearest j]
ERROR: LoadError: MethodError: no method matching merge_placeholder_expr(::Tuple{}, ::Tuple{Expr})
Closest candidates are:
  merge_placeholder_expr(::Tuple{}, ::Tuple{Vararg{Symbol,N}}) where N at /Users/roger/Documents/Workshop/QuHamiltonian/src/sub_hamiltonian.jl:5
  merge_placeholder_expr(::Tuple{Vararg{Symbol,N1}}, ::Tuple{Vararg{Symbol,N2}}) where {N1, N2} at /Users/roger/Documents/Workshop/QuHamiltonian/src/sub_hamiltonian.jl:8
Stacktrace:
 [1] macro expansion at /Users/roger/.julia/packages/MacroTools/4AjBS/src/macro.jl:74 [inlined]
 [2] QuHamiltonian.SubHamExpr(::Expr) at /Users/roger/Documents/Workshop/QuHamiltonian/src/sub_hamiltonian.jl:47
 [3] @subham(::LineNumberNode, ::Module, ::Any) at /Users/roger/Documents/Workshop/QuHamiltonian/src/sub_hamiltonian.jl:97
in expression starting at REPL[6]:1
