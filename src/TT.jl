module TT

using LinearAlgebra
using SparseArrays
using TensorOperations

import Base: size, randn

include("abstracttt.jl")
export AbstractTT,VectorTT,MatrixTT,size,randn

end # module TT
