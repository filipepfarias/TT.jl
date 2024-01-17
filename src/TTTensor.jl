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
        r1 = rounded_diagonal_index(S,e); 
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
       tmp_ps = cumsum([1;tmp_n.*tmp_r[1:tmp_d].*tmp_r[2:tmp_d+1]]);
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

function rounded_diagonal_index(S::Vector,precision::Float64)
    if (norm(S) == 0) 
        return 1
    end

    if (precision <= 0) 
        return prod(size(S))
    end

    S0 = cumsum(reverse(S) .^2 )
    ff = findall(S0 .< (precision .^ 2));
   
    return isempty(ff) ?  prod(size(S)) : prod(size(S))-ff[end]
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

"""
Tensor of all ones
   [TT]=TT_ONES(N,D), computes the d-dimensional TT-tensor of all ones 
   with mode size equal to N

   [TT]=TT_ONES(N), computes the TT-tensor of all ones with mode size 
   given in the vector N


 TT-Toolbox 2.2, 2009-2012

This is TT Toolbox, written by Ivan Oseledets et al.
Institute of Numerical Mathematics, Moscow, Russia
webpage: http://spring.inm.ras.ru/osel

For all questions, bugs and suggestions please mail
ivan.oseledets@gmail.com
---------------------------
"""
function tt_ones(n::Vector{Int})
    
    # if (numel(n) == 1)
    #     if (numel(varargin)>0)
    #         d=varargin{1}; 
    #     else
    #         d=1;
    #     end;
    #     if d>0
    #         n=n*ones(1,d);
    #     end
    # else
    #     d=numel(n);
    # end
    d = length(n);
    # if d>0
        tt=Vector{Array}(undef,d);
        for k=1:d
            tt[k] = ones(n[k]); #sqrt(n);
        end
        tt=TTTensor(tt); #Bydlocode @
    # else
    #     tt=TTTensor([[1.0]]);
    # end
    return tt
end


(-)(A::TTTensor,B::TTTensor) = A+(-1.0)*B;
(-)(B::TTTensor) = (-1.0)*B;

(+)(B::TTTensor) = B;
(+)(B::TTTensor,a::Number) = a+B;

function (+)(a::Number,B::TTTensor)
    return a*tt_ones(B.n)
end

function (+)(B::TTTensor,C::TTTensor)
    if B.d == C.d == 1
        p = getproperty.([B],propertynames(B));
        p[4] += C.core;
        return TTTensor(p...,)
    end
    
    n1=B.n;
    r1=B.r;
    r2=C.r;
    d=B.d;
    cr1=B.core;
    cr2=C.core;
    
    # Determine size of the result
    
    r = r1+r2;
    
    @assert r1[1] == r2[1] "Inconsistent sizes in mode 1."
    
    r[1]=r1[1];
    
    @assert r1[d+1] == r2[d+1] "Inconsistent sizes in mode $(d+1)"
    
    r[d+1]=r1[d+1];
    
    n=n1;
    A_d=d;
    A_n=n;
    A_r=r;
    sz=(n.*r[1:d])'*r[2:d+1];
    pos=(n.*r[1:d]) .* r[2:d+1];
    pos=cumsum([1;pos]);
    A_ps=pos;
    cr=zeros(sz);
    # if (strcmp(class(cr1),'mp'))
    # cr = mp(cr);
    # end 
    pos1=B.ps; pos2=C.ps;
    for i=1:d
        pp=zeros(r[i],n[i],r[i+1]);
        # if ( strcmp(class(cr1),'mp')) 
        #     pp = mp(pp);
        # end
        p1=cr1[pos1[i]:pos1[i+1]-1];
        p2=cr2[pos2[i]:pos2[i+1]-1];
        p1=reshape(p1,r1[i],n[i],r1[i+1]);
        p2=reshape(p2,r2[i],n[i],r2[i+1]);
        pp[1:r1[i],:,1:r1[i+1]]=p1;
        pp[r[i]-r2[i]+1:r[i],:,r[i+1]-r2[i+1]+1:r[i+1]]=p2;
        cr[pos[i]:pos[i+1]-1]=pp[:];
    end
    A_core=cr;
    
    return TTTensor(A_d,A_r,A_n,A_core,A_ps,[0])
end


function (*)(a::VecOrMat,B::TTTensor)
    d = B.d; 
    cr1 = B.core[B.ps[1]:B.ps[2]-1]; 
    cr1 = reshape(cr1, B.r[1], B.n[1]*B.r[2]);
    rnew = length(a) == 1 ? B.r[1] : size(a,1);
 
    dps = rnew*B.n[1]*B.r[2] - B.r[1]*B.n[1]*B.r[2];
    cr1 = a*cr1;
    
    # c = b;
    C_core=B.core[B.ps[2]:B.ps[d+1]-1];
    C_core = [cr1[:]; C_core[:]];
    C_ps = B.ps[:];
    C_ps[2:d+1]=C_ps[2:d+1] .+ dps;
    C_r = B.r[:];
    C_r[1] = rnew;
    
    return TTTensor(B.d,C_r,B.n,C_core,C_ps,[0])
end

(*)(a::Number,B::TTTensor) = [a]*B;

function (*)(B::TTTensor,a::VecOrMat)
    d = B.d; 
    crd = B.core[B.ps[d]:B.ps[d+1]-1]; 
    crd = reshape(crd, B.r[d]*B.n[d], B.r[d+1]);
    rnew = length(a) == 1 ? B.r[d+1] : size(a,2);
 
    ps_new = B.ps[d] + B.r[d]*B.n[d]*rnew;
    crd = crd*a;
    
    # c = b;
    C_core = B.core[:];
    if length(C_core) < ps_new-1
        append!(C_core, zeros(ps_new-1 - length(C_core)))
    end
    C_core[B.ps[d]:ps_new-1] = crd[:];
    C_core = C_core[1:ps_new-1];
    C_ps = B.ps[:];
    C_ps[d+1]= ps_new;
    C_r = B.r[:];
    C_r[d+1] = rnew;
    
    return TTTensor(B.d,C_r,B.n,C_core,C_ps,[0])
end

(*)(B::TTTensor,a::Number) = B*[a];

function norm(tt::TTTensor)
    d=tt.d;
    n=tt.n;
    r=tt.r;
    pos=tt.ps;
    cr=tt.core;
    pos1=1;
    nrm=zeros(d,1);
    core0=cr[1:r[1]*n[1]*r[2]];
    
    #Orthogonalization from left-to-tight
    for i=1:d-1
        core0=reshape(core0,r[i]*n[i],r[i+1]);
        core0,ru=qr(core0); core0 = Array(core0);
        nrm[i]=norm(ru);
        nrm[i]=max(nrm[i],1e-308);
        ru=ru./nrm[i];
        core1=cr[pos[i+1]:pos[i+2]-1];
        core1=reshape(core1,r[i+1],n[i+1]*r[i+2]);
        core1=ru*core1;
        r[i+1]=size(core0,2);
        cr[pos1:pos1-1+r[i]*n[i]*r[i+1]]=core0[:];
        cr[pos1+r[i]*n[i]*r[i+1]:pos1+r[i]*n[i]*r[i+1]+r[i+1]*n[i+1]*r[i+2]-1]=core1[:];
        core0=core1;
        pos1=pos1+r[i]*n[i]*r[i+1];
    end
    pos1=pos1+r[d]*n[d]*r[d+1]-1;
    nrm[d]=norm(core0[:]);

    return prod(nrm)
end

function dot(tt1::TTTensor,tt2::TTTensor,chunk_interval::Tuple,do_qr=false)
    
    @assert tt2.d>tt1.d "chunky dot is defined only if tt2.d>tt1.d"
    
    (chunk_start,chunk_end) = chunk_interval;

    # chunk_start and chunk_end refer to tt1
    tt2_chunk1 = chunk(tt2, chunk_start, chunk_end);
    D = dot(tt1, tt2_chunk1, do_qr);
    D = reshape(D, tt2.r(chunk_start), tt2.r(chunk_end+1));
    p = [];
    if (chunk_start>1)
        tt2_chunk1 = chunk(tt2, 1, chunk_start-1);
        p = tt2_chunk1*D;            
    end;
    if (chunk_end<tt2.d)
        tt2_chunk2 = chunk(tt2, chunk_end+1, tt2.d);
        if (isempty(p))
            p = D*tt2_chunk2;
        else
            p = tkron(p, tt2_chunk2);
        end;
    end;
    
    return p;
end

function dot(tt1::TTTensor,tt2::TTTensor,do_qr=false)
    
    if (do_qr) 
        error("Missing to implement.")
        tt1,rv1=qr(tt1, "lr");
        tt2,rv2=qr(tt2, "lr");
    end;
    
    d=tt1.d; 
    r1=tt1.r[:]; r2=tt2.r[:]; ps1=tt1.ps[:]; ps2=tt2.ps[:];
    n=tt1.n[:];
    core1=tt1.core[:]; core2=tt2.core[:];
    
    #ps is always r1(i-1)xr; but if there is a hanging thing? what to do?
    #That means, I define a dot product of two "hanging" tensors as a matrix...
    #brr.... And if it is hanging on the right? 
    # 
    # p=ones(r1(1),r2(1)); # Over r1(1) and r2(1) there is not summation blin.
    # #So, if we just sum over n(1) separatedly and leave a strange QxR1(I)xR2(I)
    # #matrix... 
    # for i=1:d
    #   cr1=core1(ps1(i):ps1(i+1)-1);
    #   cr2=core2(ps2(i):ps2(i+1)-1);
    #   p=reshape(p,[r1(i),r2(i)]);
    #   cr2=reshape(cr2,[r2(i),numel(cr2)/r2(i)]);
    #   p=p*cr2; #p is Q*r1(i)xn(i)xr2(i+1);
    #   cr1=reshape(cr1,[r1(i)*n(i),numel(cr1)/(r1(i)*n(i))]);
    #   p=reshape(p,[r1(i)*n(i),numel(p)/(r1(i)*n(i))]);
    #   p=cr1'*p;
    # end
    
    # If the first indices are not ones
    p=I(r1[1]*r2[1]);
    p = reshape(p, r1[1]*r2[1]*r1[1], r2[1]);
    
    for i=1:d
        cr1=core1[ps1[i]:ps1[i+1]-1];
        cr2=core2[ps2[i]:ps2[i+1]-1];
        cr2=reshape(cr2,r2[i], n[i]*r2[i+1]);
        
        p = p*cr2; # size r11*r21*r1-, n*r2+
        p = reshape(p,r1[1]*r2[1], r1[i]*n[i], r2[i+1]);
        p = permutedims(p, [1, 3, 2]);
        p = reshape(p, r1[1]*r2[1]*r2[i+1], r1[i]*n[i]);
        
        cr1=reshape(cr1,r1[i]*n[i], r1[i+1]);
        
        p = p*conj(cr1); # size r11*r12*r2+, r1+
        p = reshape(p, r1[1]*r2[1], r2[i+1], r1[i+1]);
        p = permutedims(p, [1, 3, 2]);
        p = reshape(p, r1[1]*r2[1]*r1[i+1], r2[i+1]);
    end;
    
    if (do_qr)
        r2old = size(rv2, 2);
        r1old = size(rv1,2);
        p = p*rv2;
        p = reshape(p, r1[1]*r2[1], r1[d+1], r2old);
        p = permutedims(p, [1, 3, 2]);
        p = reshape(p, r1[1]*r2[1]*r2old, r1[d+1]);
        p = p*conj(rv1);
        p = reshape(p, r1[1], r2[1], r2old, r1old);
        p = permutedims(p, [1,2,4,3]);
        if ( r1[1]*r2[1] == 1 ) #Save the rabbits
            p=reshape(p, r1old, r2old);
        end
    else
        p = reshape(p, r1[1], r2[1], r1[d+1], r2[d+1]);
        if ( r1[1]*r2[1] == 1 ) #Save the rabbits
            p=reshape(p, r1[d+1], r2[d+1]);
        end
    end;

    return only(p)
end
