module TT
import LinearAlgebra
import LinearAlgebra: norm, dot
import Base: size, -, +, *

include("TTTensor.jl")
export TTTensor, qtt_x, rank, ranks, rounded_diagonal_index, size, mem,+,-,*,norm,core_to_vector,vector_to_core,dot,tt_rand,tkron,+

include("TTMatrix.jl")
export  TTMatrix, vector_to_core, tt_reshape,tkron

include("solve.jl")
export amen_solve2

include("cross.jl")
export funcrs, maxvol2, rounded_diagonal_index, lu_full, reort

include("utils.jl")
export tt_ones, tt_eye, IpaS

end
