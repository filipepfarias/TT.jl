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
    tt = TTTensor(tt) #Bydlocode @
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
        e = TTMatrix([1])
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
    M = Vector{Array}(undef,d)
    M[1] = zeros(2,2,2)
    M[1][:,:,1] = [0 a;0 0]
    
    if d > 1
        M[1][:,:,2] =[0 a;0 0]
        M[d] = zeros(2,2,2)
        M[d][:,:,1] = [0 0;1 0];
        M[d][:,:,2] = [0 1;0 0];
    end
    for i=2:d-1
        M[i] = zeros(2,2,2,2);
        M[i][:,:,1,1] = I(2)
        M[i][:,:,2,1] = [0 0;1 0]
        M[i][:,:,2,2] = [0 1;0 0]
    end
    return TTMatrix(M)
end