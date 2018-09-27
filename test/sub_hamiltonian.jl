using Test, QuHamiltonian

@testset "constructor" begin
    @vertex i, j, k, l

    @subham Z[i] # => 1 instance
    @subham Z[i]Z[j] # => 2 instance
    @subham Z[i]Z[j]Z[k] # => 3 instances
    A = zeros(8, 8)
    B = zeros(4, 4)

    # multi-sites
    @subham A[i, j, k] # => 1 instance
    @subham A[i, j, k]B[j, k] # => 2 instance

    @subham(A[i, j, k]) == @subham(A[i, j, k]B[j, k])
end

@testset "vertex region" begin
    @vertex i, j

    LZ1 = @subham Z[i]
    LZ2 = @subham Z[j]
    LZ1 * LZ2 == @subham Z[i]Z[j] # i, j ∈ V1 × V2
    LZ1 * LZ1 == @subham Z[i]Z[i] # i, i ∈ V1^2
end

@testset "neighbor region" begin
    @nearest i, j, k

    LZ1 = @subham Z[i]
    LZ2 = @subham Z[j]
    LZ3 = @subham Z[k]

    @test LZ1 * LZ1 == @subham(Z[i]Z[i])
    # i, j ∈ E(V)
    @test LZ1 * LZ2 == @subham Z[i]Z[j]

    # error
    # neighbor relation should only support 2-local hamiltonians
    @test_throws ErrorException LZ1 * LZ2 * LZ3

    # this should be equivalent
    # name of place holders does not mean anything
    @test LZ1 * LZ2 == LZ2 * LZ3
end

@testset "merge placeholder" begin
    @vertex i

    A = zeros(8, 8)
    @subham(Z[i]Z[i]) == @subham($(Z*Z)[i])

    @nearest i, j, k
    @subham(A[i, j, k]A[i, j, k]) == @subham($(A * A)[i, j, k])
end
