struct TTTensor
    d::Int                       # dimension of the array
    r::Vector{Int}               # ranks of the decomposition
    n::Vector{Int}               # mode sizes of the array
    m::Vector{Int}               # mode sizes of the array
    core::Vector      # cores of the TT-decomposition store in one "long" 1D array 
    ps::Vector{Int}              # markers for position of the k-th core in array 'core'
    over::Vector{Int}
end

function TTMatrix(tt::TTTensor, m::Vector{Int}, n::Vector{Int})
    
end