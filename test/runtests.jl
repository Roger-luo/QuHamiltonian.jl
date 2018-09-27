using Test, QuHamiltonian
using Lattices

@testset "TFIM on Chain Lattice" begin

J = 0.1
h = 0.1

@nearest i, j
H = J * @sum(Z[i]Z[j]) - h * @sum(X[i])
lattice = Chain(4)
H(lattice)

end

@testset "TFIM on Square Lattice" begin

J = 0.1
h = 0.1

@nearest i, j
H = J * @sum(Z[i]Z[j]) - h * @sum(X[i])
lattice = Sqaure(4, 4)
H(lattice)

end

@testset "J1J2 on Square Lattice" begin

@neighbor 1 i, j
@neighbor 2 k, l

J1, J2 = 0.1, 0.2

H = J1 * @sum(S[i]S[j]) - J2 * @sum(S[k]S[l])
lattice = Square(4, 4)
H(lattice)

end

@testset "Heisenberg Model with rand correlation" begin

lattice = Chain(4)

@nearrest i, j
J = i, j -> rand()
H = @sum J(i, j)S[i]S[j]

lattice = Square(4, 4)
H(lattice)

end

@testset "Heisenberg Model with given correlation" begin

lattice = Chain(4)
@nearrest i, j
A = size(lattice, i, j) |> zeros
A[1, 3] = 3
J = i, j -> A[i, j]

H = @sum J(i, j)Z[i]Z[j]
H(lattice)

end

@testset "Toric Code" begin

@surround v
@face p

# A_v = @prod X[v]
# B_p = @prod Z[p]

H = -J * @sum(@prod X[v]) - J * @sum(@prod B[p])

lattice = Square(8, 8)
H(lattice)

end

@testset "SYK Model" begin

# NOTE:
# url: https://nationalmaglab.org/images/news_events/searchable_docs/winterschool/2018/Theory_2018_sachdev_syk.pdf

@region :unknow i, j, k, l

U = @parameter i, j, k, l->randn()
H = 1/(2N^(3/2)) * @sum(U[i, j, k, l]C[i]C[j]C[k]C[l]) - μ @sum(C[i] * conj(C)[j])
H(N)

end
