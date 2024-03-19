"""
% Returns a reversed TT representation
%   [TT2]=TT_REVERSE(TT,PHD)
% January 25, 2011
% Vladimir Kazeev
% vladimir.kazeev@gmail.com
% INM RAS
% Moscow, Russia
%

"""
function tt_reverse(tt,phd)

    d=size(tt,1);
    tt2=Vector{Array}(undef,d);
    
    tt2[1]=tt[d];
    for k=2:d-1
        prm=collect(1:phd+2);
        prm[phd+1]=phd+2;
        prm[phd+2]=phd+1;
        tt2[k]=permutedims(tt[d+1-k],prm);
    end
    tt2[d]=tt[1];
    
    return tt2
    end