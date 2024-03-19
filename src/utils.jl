using LinearAlgebra


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
    d = length(n)
    # if d>0
    tt = Vector{Array}(undef, d)
    for k = 1:d
        tt[k] = ones(n[k]) #sqrt(n);
    end
    tt = TTTensor(tt) #Bydlocode 
    # else
    #     tt=TTTensor([[1.0]]);
    # end
    return tt
end

function tt_eye(sz::Vector{Int64})
    d = length(sz)
    e = Vector{Array}(undef, d)
    if d != 0
        @inbounds for i in 1:d
            e[i] = I(sz[i])
        end
        e = TTMatrix(e)
    else
        e = vector_to_core(TTMatrix,[1])
    end
end

"""
A special bidiagonal matrix in the QTT-format
   [M]=IPAS(D, A)
   Generates I+a*S_{-1} matrix in the QTT-format:
   1 0 0 0
   a 1 0 0
   0 a 1 0
   0 0 a 1
 Convenient for Crank-Nicolson and time gradient matrices


 TT-Toolbox 2.2, 2009-2012

This is TT Toolbox, written by Ivan Oseledets et al.
Institute of Numerical Mathematics, Moscow, Russia
webpage: http://spring.inm.ras.ru/osel

For all questions, bugs and suggestions please mail
ivan.oseledets@gmail.com
---------------------------
"""
function IpaS(d,a)
M = Vector{Array}(undef,d);
# M[1] = zeros(2,2);
M[1]=[1 0;a 1];
if (d>1)
    # M[1](:,:,2)=[0,a;0,0];
    M[1] = cat(M[1],[0 a;0 0]; dims=3)
end;
for i=2:d-1
    M[i]=zeros(2,2,2,2);
    M[i][:,:,1,1]=I(2);
    M[i][:,:,2,1]=[0 0;1 0];
    M[i][:,:,2,2]=[0 1;0 0];
end;
if (d>1)
    M[d] = zeros(2,2,2);
    M[d][:,:,1]=I(2);
    M[d][:,:,2]=[0 0;1 0];
end;
M=TTMatrix(M); #Bydlocode
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
    d = length(n)
    res = Vector{Array{Number,3}}(undef, d)

    ni = n[1]

    for i = 2:d-1
        shabl = zeros(n[i], 2, 2)
        for j = 1:n[i]
            shabl[j, :, :] = I(2)
        end

        res[i] = shabl
        res[i][:, 2, 1] = ni .* (0:n[i]-1)
        ni = ni * n[i]
    end

    res[d] = ones(n[d], 2, 1)
    res[d][:, 2] = ni .* (0:n[d]-1)

    res[1] = ones(n[1], 2, 1)
    res[1][:, 1] = (0:n[1]-1)

    if (d == 1)
        res[1] = res[1] * [1; 0]
    end
    res = TTTensor(res) # Bydlocode ############################################
    return res
end

