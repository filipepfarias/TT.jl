module TT
import LinearAlgebra
import LinearAlgebra: norm, dot
import Base: size, -, +, *, /, reshape, getindex

include("TTTensor.jl")
export TTTensor, rank, ranks, rounded_diagonal_index, size, mem,+,-,*,norm,core_to_vector,vector_to_core,dot,tt_rand,tkron,+,getindex

include("TTMatrix.jl")
export  TTMatrix, vector_to_core, tt_reshape,tkron,*

include("solve.jl")
export amen_solve2

include("cross.jl")
export funcrs, maxvol2, rounded_diagonal_index, lu_full, reort

include("utils.jl")
export tt_ones, tt_eye, IpaS, qtt_x

include("reshape.jl")
export reshape

# Functions that must be refactored
include("to_refactor.jl")
export tt_ranks, tt_size, tt_mat_to_vec

include("core/tt_qshiftstack.jl")
export tt_qshiftstack

include("core/tt_qreshape.jl")
export tt_qreshape

include("core/tt_qtoepl.jl")
export tt_qtoepl

include("core/tt_reverse.jl")
export tt_reverse

include("core/tt_unit.jl")
export tt_unit

include("core/tt_ind2sub.jl")
export tt_ind2sub

end
