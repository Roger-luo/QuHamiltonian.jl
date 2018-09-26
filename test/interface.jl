using Revise
using Test, QuHamiltonian


@vertex i, j
@macroexpand @sum(Z[i]Z[j])
