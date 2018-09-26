struct Sum{H} end
struct LocalOperator{T, N, O, MT} end
struct SubHam{TT, OP} end
struct Vertex end

const Alias{T, MT} = Sum{SubHam{NTuple{1, Vertex}, LocalOperator{T, 1, 2, MT}}}

Alias{Float64}
