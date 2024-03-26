module TT
import LinearAlgebra
import LinearAlgebra: norm, dot, adjoint
import Base: size, -, +, *, /, reshape, getindex, setindex!,round


include("TTTensor.jl")
export TTTensor, rank, ranks, rounded_diagonal_index, size, mem,+,-,*,adjoint,norm,core_to_vector,vector_to_core,dot,tkron,getindex,setindex!,round,chunk

include("TTMatrix.jl")
export  TTMatrix, vector_to_core, tt_reshape,tkron,*,-,+,round

# include("solve.jl")
# export amen_solve2

include("solve/amen_solve.jl")
export amen_solve2, bfun3

include("solve/solve3d_2ml.jl")
export solve3d_2ml

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

include("core/tt_rand.jl")
export tt_rand

include("exp/tt_qlaplace_dd.jl")
export tt_qlaplace_dd

end
