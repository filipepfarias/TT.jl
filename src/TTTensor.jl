using LinearAlgebra
using Printf

struct TTTensor
    d::Int                       # dimension of the array
    r::Vector{Int}               # ranks of the decomposition
    n::Vector{Int}               # mode sizes of the array
    core::Vector      # cores of the TT-decomposition store in one "long" 1D array 
    ps::Vector{Int}              # markers for position of the k-th core in array 'core'
    over::Vector{Int}
end

TTTensor() = TTTensor(0,[0],[0],[0],[0],[0]);


"""
TT-tensor constructor

   T=TT_TENSOR(ARRAY) Converts from a full array with accuracy 1e-14

   T=TT_TENSOR(ARRAY,EPS) Converts from a full array with accuracy EPS

   T=TT_TENSOR(ARRAY,EPS,SZ,R1,R2) Converts from a full array which is
   treated as an array with mode sizes SZ and tail ranks R1 and R2. The
   accuracy is set to EPS

 (Adapted from) TT-Toolbox 2.2, 2009-2012 (by )

This is TT Toolbox, written by Ivan Oseledets et al.

Institute of Numerical Mathematics, Moscow, Russia

webpage: http://spring.inm.ras.ru/osel

For all questions, bugs and suggestions please mail
ivan.oseledets@gmail.com
---------------------------
"""
function TTTensor(A::Array{Number,Integer},precision=1e-14)
    n = size(A);
    d = length(n);
    r = ones(Int,d+1);
    core = zeros(0);
    c = A;
    pos=1;
    e = precision/sqrt(d-1);

    for i in 1:d-1
        m = n[i]*r[i];
        c = reshape(c,m,:);
        U,S,V = svd(c);
        # r1 = sum(cumsum(reverse(S) .^ 2) .< (e*norm(S))^2);
        # r1 = length(S) - r1;
        r1 = rounded_diagonal_index(S,e)
        U = U[:,1:r1]; S = S[1:r1]; V = V[:,1:r1];

        c = diagm(S)*V';

        r[i+1] = r1;
        append!(core,U[:])
        pos += r[i]*n[i]*r1;
    end
    append!(core,c[:]);
    
    ps=cumsum([1;n.*r[1:d].*r[2:d+1]]);

    return TTTensor(d,r,[n...],core,ps,[0])
end

function TTTensor(tt::Vector{T} where T <: AbstractArray)
    #Check if it is a canonical format
    d=prod(size(tt));
    rc=zeros(Int,d);
    all_2d=true;
    for i=1:d
       rc[i]=size(tt[i],2);
       if (size(tt[i],3) != 1 )
          all_2d = false;
       end
    end
    if ( length(unique(rc)) == 1 && all_2d)
        is_can = true;
        rc=rc[1];
    else
        is_can = false;
    end
    if ( !is_can )  
       #t=tt_tensor;
       tmp_d  = d;
       tmp_r  = ranks(tt); 
       tmp_r  = [1;tmp_r;1];
       tmp_n=[size.(tt,1)...,];  
       tmp_core=zeros(mem(tt)); 
       ps=cumsum([1;tmp_n.*tmp_r[1:d].*tmp_r[2:d+1]]);
       tmp_core[ps[1]:ps[2]-1]=tt[1];
       for i=2:d-1
          cr=tt[i]; 
          cr=permutedims(cr,[2,1,3]);
          tmp_core[ps[i]:ps[i+1]-1] = cr[:];      
       end
       cr=tt[d]; 
       cr=permutedims(cr,[2,1,3]);
       tmp_core[ps[d]:ps[d+1]-1] = cr[:];
       tmp_ps=ps;
       tmp_over=[0];
       return TTTensor(tmp_d,tmp_r,tmp_n,tmp_core,tmp_ps,tmp_over)
    else
       #t=tt_tensor;
       tmp_d=d;
       r=rc*ones(Int,d+1);
       r[1]=1; r[d+1]=1;
       tmp_r=r;
       n=zeros(Int,d);
       for i=1:d
          n[i]=size(tt[i],1);
       end
       tmp_n=n;
       tmp_ps   = cumsum([1;tmp_n.*tmp_r[1:tmp_d].*tmp_r[2:tmp_d+1]]);
       crp=zeros(tmp_ps[d+1]-1);
       cc=tt[1];
       crp[tmp_ps[1]:tmp_ps[2]-1]=cc[:];
       cc=tt[d];
       cc=transpose(cc);
       crp[tmp_ps[d]:tmp_ps[d+1]-1]=cc[:];
       for i=2:d-1
          cc=tt[i];
          cr=zeros(r[i],n[i],r[i+1]);
          for j=1:n[i]
             cr[:,j,:]=diagm(cc[j,:]);
          end
          crp[tmp_ps[i]:tmp_ps[i+1]-1]=cr[:];
       end
       tmp_core=crp;
       tmp_over=[0];
       return TTTensor(tmp_d,tmp_r,tmp_n,tmp_core,tmp_ps,tmp_over)
    end
end


"""
Computes X in the QTT-format

   [RES]=TT_X(N,D), or [RES]=TT_X(N): Vector with elements 
   0:N(1)*...*N(D)-1 is transformed into a tensor of size 
   N(1) x N(2) x ... x N(D)

   This function returns its TT-decomposition with rank 2

   This function is useful for the computation of the QTT-decomposition
   of the function y=x

   [RES]=TT_X(N,D,XMIN): Vector with elements XMIN:XMIN+N(1)*...*N(D)-1

 (Adapted from) TT-Toolbox 2.2, 2009-2012 (by )

This is TT Toolbox, written by Ivan Oseledets et al.

Institute of Numerical Mathematics, Moscow, Russia

webpage: http://spring.inm.ras.ru/osel

For all questions, bugs and suggestions please mail
ivan.oseledets@gmail.com
---------------------------
"""
function qtt_x(n::Vector{Int}) ######################## To be moved
    d = length(n);
    # if (numel(n) == 1)
    #     n=n*ones(1,d);
    # else
    #   d=numel(n);
    # end
    # xmin = 0;
    # if (numel(varargin)>1)
    #     xmin = varargin{2};
    # end;
    
    # res=cell(d,1);
    res = Vector{Array{Number,3}}(undef,d);
    
    ni=n[1];
    
    for i=2:d-1
        shabl=zeros(n[i],2,2);
        for j=1:n[i]
            shabl[j,:,:]=I(2);
        end
    
        res[i] = shabl
        res[i][:,2,1]=ni.*(0:n[i]-1);
        ni=ni*n[i];
    end
    
    res[d]=ones(n[d],2,1);
    res[d][:,2]=ni.*(0:n[d]-1);
    
    res[1]=ones(n[1],2,1);
    res[1][:,1]=(0:n[1]-1);
    
    if (d==1)
        res[1]=res[1]*[1;0];
    end;

    res = TTTensor(res); # Bydlocode ############################################
    return res
end

function rounded_diagonal_index(S::Vector,e::Float64)
    r1 = sum(cumsum(reverse(S) .^ 2) .< (e*norm(S))^2);
    return length(S) - r1;
end

function size(A::TTTensor)
    return (A.n...,)
end

size(A::TTTensor,dim::Integer) = size(A)[dim]

function ranks(cores::Vector{Array{Number,3}})
    d=size(cores,1);
    if ( d == 1 ) 
        rks=[];
        return rks
    end
    rks=zeros(Int,d-1);
    rks[1]=size(cores[1],2);

    for i=2:d-2
        rks[i] = size(cores[i],3);
    end
    rks[d-1]=size(cores[d],2);
    return rks
end

rank(A::TTTensor) = A.r;
rank(A::TTTensor,dim) = A.r[dim];

function mem(cores::Vector{Array{Number,3}})
    sz = 0;
    for cc in cores
        sz += prod(size(cc));
    end
    return sz
end

mem(A::TTTensor) = prod([ A.n .* A.r[1:A.d]; A.r[2:A.d + 1] ])