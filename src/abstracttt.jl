abstract type AbstractTT{T,N} <: AbstractArray{T,N} end

struct VectorTT{T} <: AbstractTT{T,Integer}
    cores::Vector{T}
    d::Number
    n::Vector{<:Number}
    r::Vector{<:Number}
    ps::Vector{<:Number}

    VectorTT{T}(cores::Vector{T},d::Number,n::Vector{<:Number},r::Vector{<:Number},ps::Vector{<:Number}) where {T,N} = 
    new(cores,d,n,r,ps)
end

function VectorTT(A::AbstractArray{T,N}, prec=eps()) where {T,N} 
    n = [size(A)...];
    d = prod(size(n));
    r = ones(Int,d+1);

    cores = Vector{T}();
    prec /= sqrt(d-1);
    V = [];
    for i=1:d-1
        m = n[i]*r[i];
        A = reshape(A,m,div(prod(size(A)),m));
        U,s,V = svd(A);
        r1 = trunc_sv(s,prec*norm(s));
        U = U[:,1:r1]; s = s[1:r1]; V = V[:,1:r1];
        S = Diagonal(s);
        V = V*S; A = V;
        r[i+1] = r1;
        cores = vcat(cores,U[:]);
    end
    cores = vcat(cores,V[:]);
    
    ps=cumsum([1;n.*r[1:d].*r[2:d+1]]);

    return VectorTT{T}(cores,d,n,r,ps)
end

struct MatrixTT{T} <: AbstractTT{T,Integer} 
    tt::VectorTT{T}
    m::Vector{<:Number}
    n::Vector{<:Number}

    MatrixTT{T}(tt::VectorTT{T},n::Vector{<:Number},m::Vector{<:Number}) where {T,N} = 
    new(tt,n,m)
end

function MatrixTT(A::AbstractArray{T,N}, prec=eps()) where {T,N} 
    N%2 != 0 ? throw(error("The matrix dimension must be multiple of 2, but got $N.")) : nothing;
    nm=size(A);
    d=length(nm);
    d=div(d,2);
    n=nm[1:d];
    m=nm[d+1:2*d];
    if !issparse(A)
        prm=1:2*d;
        prm=reshape(prm,d,2); 
        prm=prm';
        prm=reshape(prm,1,2*d);
        A=permutedims(A,prm);
    end;
    if (length(n.*m) == 1 )
        A=reshape(A, n.*m, 1);
    else
        A=reshape(A,n.*m);
    end
    tt=VectorTT(A,prec);
    return MatrixTT{T}(tt,[n...],[m...])
end


function trunc_sv(s,prec)
    if prec <= 0
        return length(s)
    end
    s0=cumsum(s[end:-1:1].^2);
    idx=findlast(s0 .< prec.^2);
    if isnothing(idx)
        return length(s);
    else
        return length(s)-idx;
    end
end

function size(A::AbstractTT{T,N}) where {T,N}
    sz = A.n;
    return reshape(sz, (1, length(sz)...));
end

function size(A::AbstractTT{T,N},dims) where {T,N}
    return A.n[dims];
end

function core(A::VectorTT{T},idx_core::Integer) where {T}
    iAtt = A.cores[A.ps[idx_core]:A.ps[idx_core+1]-1];
    return reshape(iAtt,[A.r[idx_core],A.n[idx_core],A.r[idx_core+1]])
end

function core(A::VectorTT{T}) where {T}
    return [core(A,idx_core) for idx_core in 1:A.d]
end

function randn(A::VectorTT{T},dir=1) where {T}
    d = copy(A.d); n = copy(A.n); r = copy(A.r); ps = copy(A.ps);
    
    if ( length(n) == 1 ) 
        n = n * ones(d);
    end
    if ( length(r) == 1 )
        r = r * ones(d+1); 
        r[1]   = 1;
        r[d+1] = 1;
    end
    
    # ps=cumsum([1;n.*r[1:d].*r[2:d+1]]);
    
    cr = zeros(T,ps[d+1]);
    if (dir>0)
        # % LR qr
        for i=1:d
            Q = randn(T,r[i]*n[i],r[i+1]);
            Q,R=qr(Q); Q = Matrix(Q);
            r[i+1] = size(Q,2);
            ps[i+1] = ps[i]+r[i]*n[i]*r[i+1];
            cr[ps[i]:(ps[i+1]-1)]=Q[:];
        end;
        cr = cr[1:(ps[d+1]-1)];
    else
        # % RL
        for i=d:-1:1
            Q = randn(T,n[i]*r[i+1], r[i]);
            Q,R=qr(Q); Q = Matrix(Q);
            Q = transpose(Q);
            r[i] = size(Q,1);
            ps[i] = ps[i+1]-r[i]*n[i]*r[i+1];
            cr[ps[i]:(ps[i+1]-1)]=Q[:];
        end;
        cr = cr[ps[1]:(ps[d+1]-1)];
        ps = ps-ps[1]+1;
    end;

    return VectorTT{T}(cr,d,n,r,ps)
end

# size(A::AbstractTT) = A.n
# numel(A::AbstractTT) = prod(size(A))