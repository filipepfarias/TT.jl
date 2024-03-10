module TT
import LinearAlgebra
import LinearAlgebra: norm, dot
import Base: size, -, +, *

include("TTTensor.jl")
export TTTensor, qtt_x, tt_ones, rank, ranks, rounded_diagonal_index, size, mem,
tt_ones,+,-,*,norm,core_to_vector, TTrand

include("solve.jl")

include("cross.jl")
export funcrs, maxvol2, rounded_diagonal_index, lu_full, reort

include("utils.jl")

end
