module QuHamiltonian

include("pauli_matrix.jl")

include("sites.jl")
include("region.jl")
include("local_operators.jl")
include("placeholder.jl")

include("sub_hamiltonian.jl")
include("symbolic.jl")
include("eval.jl")

include("conversions.jl")

# include("precompile.jl")
# _precompile_()

end # module
