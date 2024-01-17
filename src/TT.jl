module TT
import LinearAlgebra, LinearAlgebra.norm
import Base.size, Base.-, Base.+, Base.*

include("TTTensor.jl")
export TTTensor, qtt_x, tt_ones, rank, ranks, rounded_diagonal_index, size, mem,
tt_ones,+,-,*,norm

include("cross.jl")
export funcrs

include("utils.jl")

end
