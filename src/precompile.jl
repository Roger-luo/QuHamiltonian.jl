function _precompile_()
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    VecBitSites = Sites{Bit, Int, 1}
    for L = 4:10, BC in [Lattices.Periodic, Lattices.Fixed], DType in [Float32, Float64, ComplexF32, ComplexF64]
        precompile(Base.getindex, (SumSingleSiteType{DType, Matrix{DType}, Lattices.Chain{L, BC}}, VecBitSites, VecBitSites))
    end
end
