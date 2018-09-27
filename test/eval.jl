using Revise
using Test, QuHamiltonian, Lattices

# set up something useful
using LinearAlgebra, SparseArrays
⊗ = kron
II = Matrix{Float64}(I, 2, 2)
##########################

@testset "vertex sum" begin

    h = @sum Z[@vertex(i)]
    ltc = Chain(4)
    H = h(ltc)

    A = Z ⊗ II ⊗ II ⊗ II + II ⊗ Z ⊗ II ⊗ II + II ⊗ II ⊗ Z ⊗ II + II ⊗ II ⊗ II ⊗ Z

    lhs = Bit[0, 0, 0, 0]
    rhs = Bit[0, 0, 0, 0]

    for i = 1:2^4
        for j = 1:2^4
            @test H[lhs, rhs] == A[convert(Int, lhs) + 1, convert(Int, rhs) + 1]
            lhs << 1
        end
        rhs << 1
    end

    @test Matrix(H) == A
    @test sparse(H) == sparse(A)
end

@testset "edge sum" begin

    @nearest i, j
    h = @sum Z[i]Z[j]
    ltc = Chain(4)
    H = h(ltc)

    lhs = Bit[0, 0, 0, 0]
    rhs = Bit[0, 0, 0, 0]

    A = Z ⊗ Z ⊗ II ⊗ II + II ⊗ Z ⊗ Z ⊗ II + II ⊗ II ⊗ Z ⊗ Z + Z ⊗ II ⊗ II ⊗ Z
    Matrix(H)
    for i = 1:1<<4
        for j = 1:1<<4
            @test H[lhs, rhs] == A[convert(Int, lhs) + 1, convert(Int, rhs) + 1]
            lhs << 1
        end
        rhs << 1
    end

    @test Matrix(H) == A
    @test sparse(H) == sparse(A)
end

@testset "add/sub" begin

    @nearest i, j
    @vertex l
    h = @sum(X[i]X[j]) - 0.1 * @sum(Z[l])
    ltc = Chain(4)
    H = h(ltc)

    lhs = Bit[0, 0, 0, 0]
    rhs = Bit[0, 0, 0, 0]

    H[lhs, rhs]
end


# @testset "multiple vertex sum" begin
#
#     @vertex i, j
#     h = @sum Z[i]X[j] + @sum (Z * X)[i]
#     ltc = Chain(4)
#     H = h(ltc)
#
#     A =
#     # (2, 1)
#     X ⊗ Z ⊗ II ⊗ II +
#     # (3, 1)
#     X ⊗ II ⊗ Z ⊗ II +
#     # (4, 1)
#     X ⊗ II ⊗ II ⊗ Z +
#
#     # (1, 2)
#     Z ⊗ X ⊗ II ⊗ II +
#     # (3, 2)
#     II ⊗ X ⊗ Z ⊗ II +
#     # (4, 2)
#     II ⊗ X ⊗ II ⊗ Z +
#
#     # (1, 3)
#     Z ⊗ II ⊗ X ⊗ II +
#     # (2, 3)
#     II ⊗ Z ⊗ X ⊗ II +
#     # (4, 3)
#     II ⊗ II ⊗ X ⊗ Z +
#
#     # (1, 4)
#     Z ⊗ II ⊗ II ⊗ X +
#     # (2, 4)
#     II ⊗ Z ⊗ II ⊗ X +
#     # (3, 4)
#     II ⊗ II ⊗ Z ⊗ X +
#
#     # (1, 1)
#     (Z * X) ⊗ II ⊗ II ⊗ II +
#     # (2, 2)
#     II ⊗ (Z * X) ⊗ II ⊗ II +
#     # (3, 3)
#     II ⊗ II ⊗ (Z * X) ⊗ II +
#     # (4, 4)
#     II ⊗ II ⊗ II ⊗ (Z * X)
#
#     Matrix(H) == A
#
# end
