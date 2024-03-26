using LinearAlgebra

mutable struct TTMatrix
    tt::TTTensor
    n::Vector{Int}          # mode row sizes of the array
    m::Vector{Int}          # mode column sizes of the array

    function TTMatrix(tt::TTTensor, n::Vector{Int}, m::Vector{Int})
        d = tt.d
        n = (length(n) == 1) ? n .* ones(Int, d) : n
        m = (length(m) == 1) ? m .* ones(Int, d) : m

        new(tt, n, m)
    end
end

function TTMatrix(tt::TTTensor)
    n = Integer.(sqrt.(tt.n))
    m = n
    return TTMatrix(tt, n, m)
end

function TTMatrix(A::Array{T,N}, n::Vector{Int}, m::Vector{Int}, tol::Float64=1e-14) where {T<:Number,N}
    @assert (N % 2 == 0) "`A` must have a even number of dimensions."

    d = max(length(n), length(m))

    n = (length(n) < d) ? [n, ones(1, d - length(n))] : n
    m = (length(m) < d) ? [n, ones(1, d - length(m))] : m

    A = reshape(A, n..., m...)
    prm = 1:2d
    prm = reshape(prm, d, 2)
    prm = Array(prm')
    prm = reshape(prm, 1, 2d)
    permutedims!(A, A, prm)

    A = (length(n .* m) == 1) ? reshape(A, (n .* m)..., 1) : reshape(A, (n .* m)...)

    tt = TTTensor(A, tol)

    tt = !(tt.d < d) ? tt : tkron(tt, tt_ones([d - tt.d]))

    return TTMatrix(tt, n, m)
end

function TTMatrix(cc::Vector{Array})
    tt = tt_mat_to_vec(cc)
    tt = TTTensor(tt)
    d = tt.d
    n = zeros(Int, d)
    m = zeros(Int, d)
    @inbounds for i in 1:d
        n[i] = size(cc[i], 1)
        m[i] = size(cc[i], 2)
    end
    return TTMatrix(tt, n, m)
end

function TTMatrix(A::Array{T,N}, tol::Float64=1e-14) where {T<:Number,N}
    @assert (N % 2 == 0) "`A` must have a even number of dimensions."

    nm = [size(A)...,]
    d = length(nm)
    d = d ÷ 2
    n = nm[1:d]
    m = nm[d+1:2d]
    return TTMatrix(A, n, m, tol)
end

function vector_to_core(::Type{TTMatrix}, tt::Vector{T}) where {T<:AbstractArray}
    d = length(tt)
    r = zeros(Int, d + 1)
    n = zeros(Int, d)
    m = zeros(Int, d)
    ps = zeros(Int, d + 1)
    cr = Vector{Array}(undef, d)
    ps[1] = 1
    for i = 1:d
        r[i] = size(tt[i], 1)
        n[i] = size(tt[i], 2)
        m[i] = size(tt[i], 3)
        r[i+1] = size(tt[i], 4)
        cr[i] = reshape(tt[i], r[i] * n[i] * m[i] * r[i+1])
        ps[i+1] = ps[i] + r[i] * n[i] * m[i] * r[i+1]
    end
    cr = reduce(vcat, cr[:])
    tt = TTTensor(d, r, n .* m, cr, ps, [0])
    return TTMatrix(tt, n, m)
end

function core_to_vector(tm::TTMatrix)
    tt=tm.tt
    d = tt.d
    cc = Vector{AbstractArray}(undef, d)
    n = tm.n[:]
    m = tm.m[:]
    r = tt.r[:]
    ps = tt.ps[:]
    cr = tt.core[:]
    for i = 1:d
        cc[i] = reshape(cr[ps[i]:(ps[i+1]-1)], r[i], n[i], m[i], r[i+1])
    end
    return cc
end

function tkron(tt1::TTMatrix, tt2::TTMatrix)
    tt = tkron(tt1.tt, tt2.tt)
    tt3 = TTMatrix(tt, [tt1.n; tt2.n], [tt1.m; tt2.m])
    return tt3
end


function tkron(a::Vector{Any}, tt2::TTMatrix)
    @assert isempty(a) "`a` must be empty."
    return tt2
end

tkron(tt1::TTMatrix, a::Vector{Any}) = tkron(a, tt1)

function (+)(B::TTMatrix, C::TTMatrix)
    A = B.tt + C.tt
    return TTMatrix(A, B.n, B.m)
end

function (+)(a::Vector{Any}, tt2::TTMatrix)
    @assert isempty(a) "`a` must be empty."
    return tt2
end

(+)(tt1::TTMatrix, a::Vector{Any}) = (+)(a, tt1)

(-)(A::TTMatrix, B::TTMatrix) = A + (-1.0 * B);

(/)(tt1::TTMatrix, a::Number) = tt1 * (1 / a);

function (*)(A::TTMatrix, B::TTMatrix)
    d = A.tt.d
    n = A.n
    m = B.m
    p = A.m
    ra = A.tt.r
    rb = B.tt.r
    r = ra .* rb

    ps = cumsum([1; r[1:d] .* n .* m .* r[2:d+1]])
    core = Vector{Float64}(undef, ps[d+1] - 1)

    for i in 1:d
        cra = A[i]
        cra = reshape(cra, ra[i], n[i], p[i], ra[i+1])
        crb = B[i]
        crb = reshape(crb, rb[i], p[i], m[i], rb[i+1])
        cra = permutedims(cra, [1, 2, 3, 4])
        cra = reshape(cra, ra[i] * n[i] * ra[i+1], p[i])
        crb = permutedims(crb, [2, 1, 3, 4])
        crb = reshape(crb, p[i], rb[i] * m[i] * rb[i+1])
        crc::AbstractArray = cra * crb
        crc = reshape(crc, ra[i], n[i], ra[i+1], rb[i], m[i], rb[i+1])
        core[ps[i]:ps[i+1]-1] = crc[:]
    end

    C = TTTensor(d, r, n .* m, core, ps, [0])
    return TTMatrix(C, n, m)
end

function (*)(A::TTMatrix, B::TTTensor)
    n = A.n
    m = A.m
    crm = A.tt.core
    # ttm = A.tt
    psm = A.tt.ps
    d = A.tt.d
    rm = A.tt.r
    rv = B.r
    crv = B.core
    psv = B.ps
    rp = rm .* rv
    psp = cumsum([1; n .* rp[1:d] .* rp[2:d+1]])
    crp = zeros(psp[d+1] - 1)
    r = rp
    ps = psp
    for i = 1:d
        mcur = crm[psm[i]:psm[i+1]-1]
        vcur = crv[psv[i]:psv[i+1]-1]
        # make it consistent with sparse
        mcur = reshape(mcur, rm[i] * n[i] * m[i], rm[i+1])
        # mcur = mcur'
        mcur = reshape(mcur', rm[i+1] * rm[i] * n[i], m[i])
        vcur = reshape(vcur, rv[i], m[i], rv[i+1])
        vcur = permutedims(vcur, [2, 1, 3])
        vcur = reshape(vcur, m[i], rv[i] * rv[i+1])
        pcur = mcur * vcur #pcur is now rm[i+1]*rm[i]*n[i]*rv[i]*rv[i+1]
        pcur = reshape(pcur, rm[i+1], rm[i], n[i], rv[i], rv[i+1])
        #pcur=reshape(pcur,[rm[i],n[i],rm[i+1],rv[i],rv[i+1]]);
        pcur = permutedims(pcur, [2, 4, 3, 1, 5])
        #pcur=permute(pcur,[1,4,2,3,5]); 
        crp[psp[i]:psp[i+1]-1] = pcur[:]
    end
    return TTTensor(d, r, n, crp, ps, [0])
end

function (*)(A::TTMatrix, b::Array{T,N}) where {T<:Number,N}
    n = A.n
    m = A.m
    tt = A.tt
    cra = tt.core
    d = tt.d
    ps = tt.ps
    r = tt.r

    rb = size(b, 2)
    c = reshape(b, m..., rb...)

    for k = 1:d
        cr = cra[ps[k]:ps[k+1]-1]
        cr = reshape(cr, r[k], n[k], m[k], r[k+1])
        cr = permutedims(cr, [2, 4, 1, 3])
        cr = reshape(cr, n[k] * r[k+1], r[k] * m[k])
        M = prod(size(c))
        c = reshape(c, r[k] * m[k], M ÷ (r[k] * m[k]))
        c = cr * c
        M = prod(size(c))
        c = reshape(c, n[k], M ÷ n[k])
        c = permutedims(c, [2, 1])
    end
    c = c[:]
    M = prod(size(c))
    c = reshape(c, rb, M ÷ rb)
    return Array(c')
end

(*)(A::TTMatrix, b::Number) = TTMatrix(A.tt * b, A.n, A.m)
(*)(b::Number, A::TTMatrix) = TTMatrix(b * A.tt, A.n, A.m)

"""
%A=B'
    %   [A]=CTRANSPOSE(B) Compute complex conjugate transpose of a TT-matrix
    %
    %
    %
    % TT-Toolbox 2.2, 2009-2012
    %
    %This is TT Toolbox, written by Ivan Oseledets et al.
    %Institute of Numerical Mathematics, Moscow, Russia
    %webpage: http://spring.inm.ras.ru/osel
    %
    %For all questions, bugs and suggestions please mail
    %ivan.oseledets@gmail.com
    %---------------------------
"""
function adjoint(tt::TTMatrix)
    t = tt.tt
    m = tt.m
    n = tt.n  # Who wrote this code?
    d = t.d
    r = t.r
    for i = 1:d
        cr = conj(t[i]) # Blin
        cr = reshape(cr, r[i], n[i], m[i], r[i + 1])
        cr = permutedims(cr, [1, 3, 2, 4])
        t[i] = reshape(cr, r[i], m[i] * n[i], r[i + 1])
    end
    tt1 = tt
    tt1.tt = t
    tt1.m = tt.n
    tt1.n = tt.m
    return tt1
end

function core(tt1::TTMatrix, i::Integer)
    ttcr = core(tt1.tt, i)
    r = tt1.tt.r
    n = tt1.n
    m = tt1.m
    return reshape(ttcr, r[i], n[i], m[i], r[i+1])
end

"""
%Approximate TT-matrix with another one with specified accuracy
    %   [TT]=ROUND(TT,EPS,RMAX) Approximate TT-matrix with relative accuracy 
    %   EPS and maximal rank RMAX. RMAX can be array of ranks or a number
    %
    %
    %
    % TT-Toolbox 2.2, 2009-2012
    %
    %This is TT Toolbox, written by Ivan Oseledets et al.
    %Institute of Numerical Mathematics, Moscow, Russia
    %webpage: http://spring.inm.ras.ru/osel
    %
    %For all questions, bugs and suggestions please mail
    %ivan.oseledets@gmail.com
    %---------------------------
"""
round(tm::TTMatrix,eps::Float64,rmax::Vector{Int}) = round(tm.tt,eps,rmax)


"""
%Approximate TT-matrix with another one with specified accuracy
    %   [TT]=ROUND(TT,EPS) Approximate TT-matrix with relative accuracy EPS
    %
    % TT-Toolbox 2.2, 2009-2012
    %
    %This is TT Toolbox, written by Ivan Oseledets et al.
    %Institute of Numerical Mathematics, Moscow, Russia
    %webpage: http://spring.inm.ras.ru/osel
    %
    %For all questions, bugs and suggestions please mail
    %ivan.oseledets@gmail.com
    %---------------------------
"""
round(tm::TTMatrix,eps::Float64) = round(tm.tt,eps)


   
    # # [TT]=ROUND(TT,EPS)
    # # [TT]=ROUND(TT,EPS,RMAX)
    # # Approximate TT-matrix with relative accuracy EPS
    # if (nargin == 3 )
    # tt.tt=round(tt.tt,eps,rmax);
    # else
    # tt.tt=round(tt.tt,eps);    
    # end
    # return 
    # end

getindex(tt::TTMatrix, i::Integer) = core(tt, i)