using Test, QuHamiltonian

@testset "neighbor" begin
    @neighbor 2 i,j
    @test i == QuHamiltonian.Neighbor{2}(:i)
    @test j == QuHamiltonian.Neighbor{2}(:j)
end

@testset "nearest" begin
    @nearest i, j
    @test i == QuHamiltonian.Neighbor{1}(:i)
    @test j == QuHamiltonian.Neighbor{1}(:j)
end

@testset "vertex" begin
    @vertex i, j
    @test i == QuHamiltonian.Vertex(:i)
    @test j == QuHamiltonian.Vertex(:j)
end

@testset "surround" begin
    @surround i, j
    @test i == QuHamiltonian.Surround(:i)
    @test j == QuHamiltonian.Surround(:j)

    # NOTE: surround a vertex place holder is more readable?
    # @vertex k
    # @surround k
end

@testset "face" begin
    @face i, j
    @test i == QuHamiltonian.Face(:i)
    @test j == QuHamiltonian.Face(:j)
end
