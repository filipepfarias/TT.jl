module TT
import LinearAlgebra
import LinearAlgebra: norm, dot
import Base: size, -, +, *

include("TTTensor.jl")
export TTTensor, qtt_x, tt_ones, rank, ranks, rounded_diagonal_index, size, mem,
tt_ones,+,-,*,norm

include("cross.jl")
export funcrs

include("utils.jl")

end
