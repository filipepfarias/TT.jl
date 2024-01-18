"""
Cross approximation of a function of a TT-tensor, Method 1
    
    [Y]=FUNCRS(TT,FUN,EPS,Y,NSWP)
    
Computes approximation to the function FUN(TT) with accuracy EPS
Auxiliary parameters:  Y (initial approximation), NSWP 
(number of sweeps in the method 
Much faster then usual cross by vectorized computation of subtensors
(Adapted from) TT-Toolbox 2.2, 2009-2012 (by )

This is TT Toolbox, written by Ivan Oseledets et al.
Institute of Numerical Mathematics, Moscow, Russia
webpage: http://spring.inm.ras.ru/osel

For all questions, bugs and suggestions please mail
ivan.oseledets@gmail.com
---------------------------
"""
function funcrs(tt,fun,precision,y,n_sweeps)
    dropsweeps = 4;
    kick_rank=10;
    verb=true;
    # if (~isempty(y))
    #    yold=y;
    # end
    #For this procedure we need only local (!) indices, since it
    #is more or less equivalent to orthogonalization; 
    d=tt.d; 
    ps=tt.ps[:];
    core=tt.core[:];
    n=tt.n[:];
    r=tt.r[:];
    
    ry=y.r[:];
    psy=y.ps[:];
    cry=y.core[:];
    phi= Vector{Any}(undef,d+1);
    phi[d+1]=1;
    phi[1]=1;
    
    #Warmup procedure: orthogonalize from right-to-left & maxvol
    pos1=psy[d];
    #Cores i
    for i in d-1:-1:1
        cr=cry[pos1:pos1+ry[i+1]*n[i+1]*ry[i+2]-1];
        cr2=cry[pos1-ry[i]*n[i]*ry[i+1]:pos1-1];
        cr2=reshape(cr2,ry[i]*n[i],ry[i+1]);
        cr=reshape(cr,ry[i+1],n[i+1]*ry[i+2]);;
        cr=cr';
        cr,rm=qr(cr); cr = Array(cr); 
        ry[i+1]=size(cr,2);  
        ind=maxvol2(cr); 
        r1=cr[ind,:];
        cr=cr/r1; 
        cr=cr';  
        cry[pos1:pos1+ry[i+1]*n[i+1]*ry[i+2]-1]=cr[:];
        pos1=pos1-ry[i]*n[i]*ry[i+1];
        cr2=cr2*(r1*rm)';  
        cry[pos1:pos1+ry[i]*n[i]*ry[i+1]-1]=cr2[:]; 
        #Take phi matrix; convolve from right with current cors of V
        cr0=core[ps[i+1]:ps[i+2]-1];
        cr0=reshape(cr0,r[i+1]*n[i+1],r[i+2]); 
        cr0=cr0*phi[i+2]; #cr0 is now r(i)*n(i)*ry(i+1);
        cr0=reshape(cr0,r[i+1],n[i+1]*ry[i+2]);
        phi[i+1]=cr0[:,ind]; 
        
        #  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
        #tmpy_core=cry[:];
        #tmpy_r=ry[:];
        #tmpy_ps=psy[:];
        #keyboard
        
        #Orthogonalization & maxvol is "local operation" (touches two
        #cores)
        #Computation of elements requries storage of phi matrices   
    end
    
    #Truncate cry
    #pos1=pos1-n(1)*r(2);
    cry=cry[pos1:end]; 
    tmpy_core=cry[:];
    tmpy_r=ry[:];
    tmpy_ps=psy[:];
    #keyboard
    swp=1;
    
    
    yold=[];
    not_converged = true;
    pos1=1;
    last_sweep = false;
    while ( swp < n_sweeps && not_converged )

        max_er=0;
        cry_old=cry[:]; 
        psy=cumsum([1;n.*ry[1:d].*ry[2:d+1]]); 
        pos1=1;
        
        for i=1:d-1
            #fprintf('i=%d \n',i);
            #We care for two cores, with number i & number i+1, and use
            #psi(i) and psi(i+2) as a basis; also we will need to recompute
            #psi(i+1)
            ps1=phi[i]; ps2=phi[i+2];
            #Take current core; convolve it with 
            #core=core(ps
            #Compute (!) superblock and function (!) of it
            
            cr1=core[ps[i]:ps[i+1]-1];
            cr2=core[ps[i+1]:ps[i+2]-1];
            cr1=reshape(cr1,r[i],n[i]*r[i+1]); 
            cr1=ps1*cr1;
            cr2=reshape(cr2,r[i+1]*n[i+1],r[i+2]);
            cr2=cr2*ps2;
            cr1=reshape(cr1,ry[i]*n[i],r[i+1]);
            cr2=reshape(cr2,r[i+1],n[i+1]*ry[i+2]);
            cr=cr1*cr2;
            cr=reshape(cr,ry[i]*n[i],n[i+1]*ry[i+2]);
            cr=fun.(cr);  #Elements are evaluated here!
            cr=reshape(cr,ry[i]*n[i],n[i+1]*ry[i+2]);
            #Check for local approximation of cr for the error
            cry1=cry[pos1:pos1+ry[i]*n[i]*ry[i+1]-1]; 
            
            cry2=cry_old[psy[i+1]:psy[i+2]-1]; 
            
            cry1=reshape(cry1,ry[i]*n[i],ry[i+1]);
            cry2=reshape(cry2,ry[i+1],n[i+1]*ry[i+2]);
            appr=cry1*cry2; 
            er=norm(appr-cr)/norm(cr);
            max_er=max(er,max_er);
            #Compute SVD of cr
            
            if ((swp % dropsweeps)!=0 && (swp > 1) && !last_sweep)
                u,s,v=svd(cr-appr); 
            else
                u,s,v=svd(cr); 
            end;
            # s=diag(s);
            r2=rounded_diagonal_index(s,precision*norm(cr)/sqrt(d-1)); 
            s=s[1:r2]; u=u[:,1:r2]; v=v[:,1:r2]; 
            
            v=conj(v)*diagm(s);
            
            if ((swp % dropsweeps)!=0 && (swp > 1) && !last_sweep)
                u = [cry1 u]; 
                v = [cry2' v];
                u,rv=qr(u); u = Array(u); 
                v = v*rv';         
            else
                if !last_sweep
                    #Kick rank of u
                    ur=randn(size(u,1),kick_rank); 
                    #Orthogonalize ur to u by Golub-Kahan reorth
                    u=reort(u,ur); 
                    radd=size(u,2)-r2;
                    if radd > 0 
                        vr=zeros(size(v,1),radd);
                        v=[v vr];
                    end
                end;
            end;
            
            #vr=randn(size(v,1),kick_rank);   
            #v=reort(v,vr);
            #radd=size(v,2)-rnew; 
            #if ( radd > 0 )
            #    ur=zeros(size(u,1),radd);
            #    u=[u,ur];
            #end
            #rnew=rnew+radd;
            
            
            
            
            #u=[u,uadd]; v=[v,vadd];
            #[u,ru]=qr(u,0); 
            #v=v*(ru.');
            r2=size(u,2);
            
            ind=maxvol2(u); 
            r1=u[ind,:]; 
            u=u/r1; v=v*r1'; v=v';
            #Compute new phi
            ry[i+1]=r2;
            if length(cry) < pos1+ry[i]*n[i]*ry[i+1]-1 ## Substitutes MATLAB's bad allocation
                append!(cry,zeros(pos1+ry[i]*n[i]*ry[i+1]-1 - length(cry)))
            end
            
            cry[pos1:pos1+ry[i]*n[i]*ry[i+1]-1]=u[:];
            pos1=pos1+ry[i]*n[i]*ry[i+1];
            if length(cry) < pos1+ry[i+1]*n[i+1]*ry[i+2]-1 ## Substitutes MATLAB's bad allocation
                append!(cry,zeros(pos1+ry[i+1]*n[i+1]*ry[i+2]-1 - length(cry)))
            end
            cry[pos1:pos1+ry[i+1]*n[i+1]*ry[i+2]-1]=v[:];
            
            #Now we have to: cr1 with phi from the left 
            phi[i+1]=cr1[ind,:]; #phi[i+1]=phi[i+1]; 
        end; 
        
        #Truncate local memory
        cry=cry[1:pos1+ry[d]*n[d]*ry[d+1]-1];
        psy=cumsum([1;n.*ry[1:d].*ry[2:d+1]]);
        tmpy_core=cry[:];
        tmpy_r=ry[:];
        tmpy_ps=psy[:];
        #keyboard;
        cry_old=cry[:];
        #m=numel(cry);
        #cry=[zeros(size(cry)),cry];
        #pos1=pos1+m;
        #cry2=cry_old(psy(i+1):psy(i+2)-1);
        cry=cry[psy[d]:psy[d+1]-1];  #Start--only two cores
        for i=d-1:-1:1
            #fprintf('i=%d \n',i);
            #We care for two cores, with number i & number i+1, and use
            #psi(i) and psi(i+2) as a basis; also we will need to recompute
            #psi(i+1)
            ps1=phi[i]; ps2=phi[i+2];
            #Take current core; convolve it with 
            #core=core(ps
            #Compute (!) superblock and function (!) of it
            cr1=core[ps[i]:ps[i+1]-1];
            cr2=core[ps[i+1]:ps[i+2]-1];
            cr1=reshape(cr1,r[i],n[i]*r[i+1]);
            cr1=ps1*cr1;
            cr2=reshape(cr2,r[i+1]*n[i+1],r[i+2]);
            cr2=cr2*ps2;
            cr1=reshape(cr1,ry[i]*n[i],r[i+1]);
            cr2=reshape(cr2,r[i+1],n[i+1]*ry[i+2]);
            cr=cr1*cr2;
            cr=reshape(cr,ry[i]*n[i],n[i+1]*ry[i+2]);
            cr=fun.(cr); #Elements are evaluated here!
            cr=reshape(cr,ry[i]*n[i],n[i+1]*ry[i+2]);
            #Check for local approximation of cr for the error
            #cry1=cry(pos1-ry(i)*n(i)*ry(i+1):pos1-1);
            cry1=cry_old[psy[i]:psy[i+1]-1];
            #cry2=cry(ry(:ry(i+1)*n(i+1)*ry(i+2));
            #cry1=cry(1:ry(i)*n(i)*ry(i+1));
            #cry2=cry(ry(i)*n(i)*ry(i+1)+1:ry(i)*n(i)*ry(i+1)+ry(i+1)*n(i+1)*ry(i+2));
            cry2=cry[1:ry[i+1]*n[i+1]*ry[i+2]];
            cry = cry[ry[i+1]*n[i+1]*ry[i+2]+1:end];
            #cry2=cry_old(psy(i+1):psy(i+2)-1);
            #cry(1:(ry(i+1)*n(i+1)*ry(i+2)))=[]; #Delete both of first cores
            #cry(1:ry(i)*n(i)*ry(i+1)+ry(i+1)*n(i+1)*ry(i+2))=[];
            cry1=reshape(cry1,ry[i]*n[i],ry[i+1]);
            cry2=reshape(cry2,ry[i+1],n[i+1]*ry[i+2]);
            appr=cry1*cry2;
            er=norm(appr-cr)/norm(cr);
            #er
            max_er=max(er,max_er); 
            
            if (mod(swp,dropsweeps)!=0)&&(swp>1)&&(!last_sweep)
                u,s,v=svd(cr-appr);
            else
                #Compute SVD of cr
                u,s,v=svd(cr);
            end;
            #s=diag(s);
            r2=rounded_diagonal_index(s,precision*norm(cr)/sqrt(d-1));
            s=s[1:r2]; u=u[:,1:r2]; v=conj(v[:,1:r2]);
            #Make it standard
            u=u*diagm(s); 
            
            
            if (mod(swp,dropsweeps)!=0)&&(swp>1)&&(!last_sweep)
                u = [cry1 u]; 
                v = [cry2' v];
                v,rv=qr(v); v = Array(v);
                u = u*rv'; 
            else
                if (!last_sweep)
                    #Kick rank
                    vr=randn(size(v,1),kick_rank);
                    v=reort(v,vr);
                    radd=size(v,2)-r2;
                    if ( radd > 0 )
                        ur=zeros(size(u,1),radd);
                        u=[u ur];  
                    end
                end;
            end;
            
            
            
            r2=size(v,2);
            ind=maxvol2(v);
            r1=v[ind,:]; 
            v=v/r1; v=v';
            u=u*r1';
            #Compute new phi;
            ry[i+1]=r2;
            u=u[:]; u=u'; v=v[:]; v=v'; 
            # cry=[u v cry]; 
            prepend!(cry,u,v);
            #keyboard;
            #We need new memory for
            #cry(pos1:pos1+ry(i+1)*n(i+1)*ry(i+2)-1)=v(:); #Here memory has to be (?) enlarged 
            #if ( pos1 <= ry(i)*n(i)*ry(i+1) ) #Have to enlarge memory
            #  cry=cry(pos1:numel(cry));
            #  cry=[u(:)',cry];
            #  pos1=ry(i)*n(i)*ry(i+1)+1;
            #else
            #   pos1=pos1-ry(i)*n(i)*ry(i+1);
            #  cry(pos1:pos1+ry(i)*n(i)*ry(i+1)-1)=u(:);
            #end
            #Now we have to: cr1 with phi from the left 
            phi[i+1]=cr2[:,ind]; phi[i+1]=phi[i+1];
        end;
        #Truncate local memory
        #cry=cry(pos1:numel(cry));
        
        # keyboard;
        psy=cumsum([1;n.*ry[1:d].*ry[2:d+1]]);
        
        tmpy_core=cry[:]; 
        tmpy_r=ry[:];
        tmpy_ps=psy[:];
        #keyboard;
        # if ( isempty(yold) )
        if ( !isa(yold,TTTensor) )
            yold=y;
            er_nrm=1;
        else
            
            er_nrm=norm(yold-y)/norm(y);
            yold=y;
        end
        if ( verb )
            msg = @sprintf "sweep=%d, er=%3.2e er_nrm=%3.2e \n" swp max_er er_nrm;
            print(msg)
        end
        if (last_sweep)
            break;
        end;
        if (er_nrm<precision)
            last_sweep=true;
        end; 
        swp+=1;
    end     
    psy=cumsum([1;n.*ry[1:d].*ry[2:d+1]]);
    
    tmpy_core=cry[:];
    tmpy_r=ry[:];
    tmpy_ps=psy[:];
    
    return TTTensor(y.d,tmpy_r,y.n,tmpy_core,tmpy_ps,y.over)
end

"""
Maximal volume submatrix in an tall matrix
[IND]=MAXVOL2(A,IND) Computes maximal volume submatrix in A starting from IND
[IND]=MAXVOL2(A) Computes maximal volume submatrix in A starting from LU
Returns rows indices that contain maximal volume submatrix

Reference: 
How to find a good submatrix / S.A. Goreinov, 
I.V. Oseledets, D.V. Savostyanov, E.E. Tyrtyshnikov, N.L. Zamarashkin
// Matrix Methods: Theory, Algorithms, Applications / 
Ed. by V. Olshevsky, E. Tyrtyshnikov. ? World Scientific 
Publishing, 2010. ? Pp. 247?256.


(Adapted from) TT-Toolbox 2.2, 2009-2012 (by )

This is TT Toolbox, written by Ivan Oseledets et al.
Institute of Numerical Mathematics, Moscow, Russia
webpage: http://spring.inm.ras.ru/osel

For all questions, bugs and suggestions please mail
ivan.oseledets@gmail.com
--------------------------- 
"""
function maxvol2(a,do_qr=false,do_lu_full=false,niters=100,precision=5e-2)
    
    n=size(a,1); r=size(a,2);
    if ( n <= r )
        ind = 1:n;
        ## Not tested ##
        return ind
    end
    
    if ( do_qr ) 
        a,_ = qr(a); a = Array(a); 
    end
    
    
    #Initialize
    if ( do_lu_full )
        ind = lu_full(a); println("lu_full");
    else
        _,_,p=lu(a);
        ind = p[1:r];
    end    
    sbm = a[ind,:];
    b = a / sbm;
    
    #Start iterations
    iter=0;
    while (iter <= niters);
        #[mx0,big_ind] = max(abs(b(:)));
        mx0, big_ind = findmax( abs.(b) );
        #[i0,j0] = ind2sub([n,r],big_ind);
        if ( mx0 <= 1 + precision ) 
            ind=sort(ind);
            return ind
        end 
        k = ind[big_ind[2]];  #This is the current row that we are using
        b = b + b[:,big_ind[2]]*(b[k,:]'-b[big_ind[1],:]')/b[big_ind]; 
        ind[big_ind[2]] = big_ind[1];
    end

    return ind
end

"""
LU with full pivoting
[U,V,IND] = LU_FULL(A)

Reference: 
How to find a good submatrix / S.A. Goreinov, 
I.V. Oseledets, D.V. Savostyanov, E.E. Tyrtyshnikov, N.L. Zamarashkin
// Matrix Methods: Theory, Algorithms, Applications / 
Ed. by V. Olshevsky, E. Tyrtyshnikov. ? World Scientific 
Publishing, 2010. ? Pp. 247?256.


TT-Toolbox 2.2, 2009-2012

(Adapted from) TT-Toolbox 2.2, 2009-2012 (by )

This is TT Toolbox, written by Ivan Oseledets et al.
Institute of Numerical Mathematics, Moscow, Russia
webpage: http://spring.inm.ras.ru/osel

For all questions, bugs and suggestions please mail
ivan.oseledets@gmail.com
--------------------------- 
"""
function lu_full(a)
    b = a;
    n, r = size(b);
    ind = zeros(r);
    
    for i in 1:r 
        mx0, big_ind = findmax( abs.(b) );
        b = b - b[:,big_ind[2]]*b[big_ind[1],:]/b[big_ind];
        ind[big_ind[2]] = i0;
    end
    return ind
end

"""
Golub-Kahan reorthogonalization
[Y]=REORT(U,UADD) Given orthonormal matrix U, orthonormalize UADD to it
to get orthogonal basis in [U,UADD]. Faster then the QR-decomposition,
using the Golub-Kahan method


(Adapted from) TT-Toolbox 2.2, 2009-2012 (by )

This is TT Toolbox, written by Ivan Oseledets et al.
Institute of Numerical Mathematics, Moscow, Russia
webpage: http://spring.inm.ras.ru/osel

For all questions, bugs and suggestions please mail
ivan.oseledets@gmail.com
--------------------------- 
"""
function reort(u,uadd)
    usave=u;
    uadd_save=uadd;
    if (size(uadd,2)==0)
        y = u;
        return y;
    end;
    
    if (size(u,1) == size(u,2) )
        y=u;
        return y
    end
    
    if (size(u,2) + size(uadd,2) >= size(u,1) )
        uadd=uadd[:,1:(size(u,1)-size(u,2))];
    end
    
    mvr=u'uadd; unew=uadd-u*mvr; 
    reort_flag=true;
    j=1;
    # z=[];
    while (reort_flag && j <= 20 )
        # Nice vectorization!
        norm_unew=sum(unew.^2,dims=1);
        norm_uadd=sum(uadd.^2,dims=1);
        # if (sum(norm_unew,2)<1e-28*sum(norm_uadd,2))
        #     unew = [];
        #     break;
        # end;
        # reort_flag=~isempty(find(norm_unew <= 0.25*norm_uadd,1));
        reort_flag=any(norm_unew .<= 0.25*norm_uadd);
        
        unew,_=qr(unew); unew = Array(unew) #Here it is ok. TO BE verified ###################################
        if (reort_flag)
            uadd=unew;
            su=u'unew;
            unew=unew-u*su;
            #Kill machine zeros
            #unew(unew==0)=max(1e-308,norm(unew,'fro')*1e-15);
            j=j+1;
            #if ( j == 3 && isempty(z) ) #Kill machine zeros
            #   z=randn(size(unew,1),1); z=z/norm(z);
            #   u = u - 2*z*(z'*u);
            #   unew = unew - 2*z*(z'*unew);
            #end
        end
    end
    if ( reort_flag )
        display("Reort failed to reort!");
        y,_ =qr([u unew]); y = Array(y)
    else
        y=[u unew];
    end
    #if ( ~isempty(z) ) 
    #   y = y - 2*z*(z'*y);
    #end
    return y
end

function [y]=multifuncrs(X, funs, eps, varargin)
    # Cross approximation of a (vector-)function of several TT-tensors.
    #   [Y]=MULTIFUNCRS(X,FUNS,EPS, VARARGIN)
    #   Computes approximation to the functions FUNS(X{1},...,X{N}) with accuracy EPS
    #   X should be a cell array of nx TT-tensors of equal sizes.
    #   The function FUNS should receive a 2d array V of sizes I x N, where the
    #   first dimension stays for the reduced set of spatial indices, and  the
    #   second is the enumerator of X.
    #   The returned sizes should be I x D2, where D2 is the number of
    #   components in FUNS. D2 should be either provided as the last (d+1)-th
    #   TT-rank of the initial guess, or given explicitly as an option (see
    #   below).
    #   For example, a linear combination reads FUNS=@(x)(x*W), W is a N x D2
    #   matrix.
    #
    #   Options are provided in form
    #   'PropertyName1',PropertyValue1,'PropertyName2',PropertyValue2 and so
    #   on. The parameters are set to default (in brackets in the following)
    #   The list of option names and default values are:
    #       o y0 - initial approximation [random rank-2 tensor]
    #       o nswp - maximal number of DMRG sweeps [10]
    #       o rmax - maximal TT rank [Inf]
    #       o verb - verbosity level, 0-silent, 1-sweep info, 2-block info [1]
    #       o kickrank - the rank-increasing parameter [5]
    #       o d2 - the last rank of y, that is dim(FUNS) [1]
    #       o qr - do (or not) qr before maxvol [false]
    #
    #   The method is based on the alternating approximation, with 
    #   the one-block enrichment via KICKRANK random vectors or randomized AMR.
    #   
    #
    # TT-Toolbox 2.2, 2009-2012
    #
    #This is TT Toolbox, written by Ivan Oseledets et al.
    #Institute of Numerical Mathematics, Moscow, Russia
    #webpage: http://spring.inm.ras.ru/osel
    #
    #For all questions, bugs and suggestions please mail
    #ivan.oseledets@gmail.com
    #---------------------------
    
    nswp = 10;
    kickrank = 5;
    y = [];
    verb = 1;
    # kicktype = 'rand';
    kicktype = "amr-two";
    pcatype = "svd";
    #pcatype = "uchol";
    rmax = Inf;
    d2 = 1;
    wasrand = false;
    trunctype = "fro";
    do_qr = false;
    # trunctype = 'cross';
    
    # for i=1:2:length(varargin)-1
    #     switch lower(varargin{i})
    #         case 'nswp'
    #             nswp=varargin{i+1};
    #         case 'y0'
    #             y=varargin{i+1};
    #         case 'kickrank'
    #             kickrank=varargin{i+1};
    #         case 'rmax'
    #             rmax=varargin{i+1};            
    #         case 'verb'
    #             verb=varargin{i+1};
    #         case 'kicktype'
    #             kicktype=varargin{i+1};            
    #         case 'pcatype'
    #             pcatype=varargin{i+1};
    #         case 'trunctype'
    #             trunctype=varargin{i+1};            
    #         case 'd2'
    #             d2=varargin{i+1};            
    #         case 'qr'
    #             do_qr = varargin{i+1};
    #         otherwise
    #             error('Unrecognized option: %s\n',varargin{i});
    #     end
    # end
    
    nx = numel(X);
    d = X{1}.d;
    n = X{1}.n;
    rx = zeros(d+1,nx);
    crX = cell(d,nx);
    for i=1:nx
        rx(:,i) = X{i}.r;
        crX(:,i) = core2cell(X{i});
    end;
    
    if (isempty(y))
        ry = d2*ones(d+1,1); ry(1)=1;
        y = tt_rand(n, d, ry);
        wasrand = true;
    end;
    ry = y.r;
    cry = core2cell(y);
    
    Ry = cell(d+1,1);
    Ry{1} = 1; Ry{d+1}=1;
    Rx = cell(d+1,nx);
    for i=1:nx
        Rx{1,i}=1; Rx{d+1,i}=1;
    end;
    
    block_order = [+(d), -(d)];
    
    # Orth
    for i=1:d-1
        cr = cry{i}; # r1,n,d2,r2
        cr = reshape(cr, ry(i)*n(i), ry(i+1));
        [cr, rv]=qr(cr, 0);    
        cr2 = cry{i+1};
        cr2 = reshape(cr2, ry(i+1), n(i+1)*ry(i+2));
        cr2 = rv*cr2;
        ry(i+1) = size(cr, 2);
        cr = reshape(cr, ry(i), n(i), ry(i+1));
        cry{i+1} = reshape(cr2, ry(i+1), n(i+1), ry(i+2));
        cry{i} = cr;
    
        # Interface matrix for Y        
        Ry{i+1} = Ry{i}*reshape(cr, ry(i), n(i)*ry(i+1));
        Ry{i+1} = reshape(Ry{i+1}, ry(i)*n(i), ry(i+1));
        if (wasrand)
            curind = [];
            while numel(curind)<ry(i+1)
                curind = [curind; ceil(rand(ry(i+1), 1)*(n(i)*ry(i)))];
                curind = unique(curind);
            end;
            curind = curind(1:ry(i+1));
        else
            
            curind = maxvol2(Ry{i+1},"qr",do_qr);
        end;    
        Ry{i+1} = Ry{i+1}(curind, :);
        # Interface matrices for X
        for j=1:nx
            Rx{i+1,j} = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
            Rx{i+1,j} = Rx{i,j}*Rx{i+1,j};
            Rx{i+1,j} = reshape(Rx{i+1,j}, ry(i)*n(i), rx(i+1,j));
            Rx{i+1,j} = Rx{i+1,j}(curind, :);
        end;   
    end;
    
    
    d2 = ry(d+1);
    ry(d+1) = 1;
    cry{d} = permute(cry{d}, [3,1,2]); # d2, rd, nd
    
    last_sweep = false;
    swp = 1;
    
    # dy_old = ones(d,1);
    dy = zeros(d,1);
    max_dy = 0;
    # For extra-rank addition
    # dpows = ones(d,1)*min_dpow;
    # dranks = zeros(d,1);
    
    cur_order = block_order;
    order_index = 2;
    i = d;
    dir = sign(cur_order(order_index));
    
    # DMRG sweeps
    while (swp<=nswp)||(dir>0)
        
        oldy = reshape(cry{i}, d2*ry(i)*n(i)*ry(i+1), 1);
        
        if (~last_sweep)
            # Compute the X superblocks
            curbl = zeros(ry(i)*n(i)*ry(i+1), nx);
            for j=1:nx
                cr = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                cr = Rx{i,j}*cr;
                cr = reshape(cr, ry(i)*n(i), rx(i+1,j));
                cr = cr*Rx{i+1,j};
                curbl(:,j) = cr(:);
            end;
            # Call the function
            newy = funs(curbl); # sizes: rnr x nx -> rnr x d2
            # Multiply with inverted Ry
            newy = reshape(newy, ry(i), n(i)*ry(i+1)*d2);
            newy = (Ry{i}) \ newy;
            newy = reshape(newy, ry(i)*n(i)*ry(i+1), d2);
            newy = reshape(newy.', d2*ry(i)*n(i), ry(i+1));
            newy = newy / (Ry{i+1});
            newy = reshape(newy, d2*ry(i)*n(i)*ry(i+1), 1);
        else
            newy = oldy;
        end;
    
        dy(i) = norm(newy(:)-oldy)/norm(newy(:));
        max_dy = max(max_dy, dy(i));
    
        # Truncation
        if (dir>0) # left-to-right
            newy = reshape(newy, d2, ry(i)*n(i)*ry(i+1));
            newy = reshape(newy.', ry(i)*n(i), ry(i+1)*d2);
        else
            newy = reshape(newy, d2*ry(i), n(i)*ry(i+1));
        end;
    
        if (kickrank>=0)
            [u,s,v]=svd(newy, "econ");
            s = diag(s);
            if (strcmp(trunctype, "fro"))||(last_sweep)
                r = my_chop2(s, eps/sqrt(d)*norm(s));
            else            
                # Truncate taking into account the (r+1) overhead in the cross
                cums = (s.*(2:numel(s)+1)').^2;
                cums = cumsum(cums(end:-1:1));
                cums = cums(end:-1:1)./cums(1);
                r = find(cums<(eps^2/d), 1);
                if (isempty(r))
                    r = numel(s);
                end;
            end;
            r = min(r, rmax);
            r = min(r, numel(s));
        else
            if (dir>0)
                [u,v]=qr(newy, 0);
                v=v';
                r = size(u,2);
                s = ones(r,1);
            else
                [v,u]=qr(newy.', 0);
                v=conj(v);
                u=u.';
                r = size(u,2);
                s = ones(r,1);
            end;
        end;
    
        if (verb>1)
            fprintf("=multifuncrs=   block %d{%d}, dy: %3.3e, r: %d\n", i, dir, dy(i), r);
        end;    
        
        # Kicks and interfaces
        if (dir>0)&&(i<d) # left-to-right, kickrank, etc
            u = u(:,1:r);
            v = v(:,1:r)*diag(s(1:r));
    
            # kick
            radd = 0; rv = 1;
            if (~last_sweep)&&(kickrank>0)
                if (strcmp(kicktype, "amr-two"))
                    # AMR(two)-like kick. See also the M.Sc.Thesis by D. Zheltkov.
                    # The left indices are nested, but the right are chosen
                    # randomly. In Zheltkov's work, from all possible n^(d-k)
                    # values. However, in the functional-cross it would result
                    # in a d^2 complexity. Here, I use only the neighbouring
                    # core for randomization. Additionally, the actual kick is
                    # performed via the Z=PCA(supercore), since otherwise the
                    # rank grows too high.
                    
                    # Compute the X superblocks
                    ind2 = unique(ceil(rand(ry(i+1), 1)*(ry(i+2)*n(i+1))));
                    rkick = numel(ind2);
                    curbl = zeros(ry(i)*n(i)*rkick, nx);
                    for j=1:nx
                        cr1 = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                        cr1 = Rx{i,j}*cr1;
                        cr1 = reshape(cr1, ry(i)*n(i), rx(i+1,j));
                        cr2 = reshape(crX{i+1,j}, rx(i+1,j)*n(i+1), rx(i+2,j));                    
                        cr2 = cr2*Rx{i+2,j}; # now its size rx
                        cr2 = reshape(cr2, rx(i+1,j), n(i+1)*ry(i+2));
                        cr2 = cr2(:, ind2);
                        curbl(:,j) = reshape(cr1*cr2, ry(i)*n(i)*rkick, 1);
                    end;
                    # Call the function
                    uk = funs(curbl); # rnr, d2
                    uk = reshape(uk, ry(i), n(i)*rkick*d2);
                    uk = Ry{i} \ uk;
                    uk = reshape(uk, ry(i)*n(i), rkick*d2);
                    if (strcmp(pcatype, "svd"))
                        [uk,sk,vk]=svd(uk, "econ");
                        uk = uk(:,1:min(kickrank, size(uk,2)));
                    else
                        uk = uchol(uk.', kickrank+1);
                        uk = uk(:,end:-1:max(end-kickrank+1,1));
                    end;
                else
                    uk = rand(ry(i)*n(i), kickrank);
                end;
                [u,rv]=qr([u,uk], 0);
                radd = size(uk,2);
            end;
            v = [v, zeros(ry(i+1)*d2, radd)];
            v = rv*(v');
            r = size(u,2);
    
            cr2 = cry{i+1};
            cr2 = reshape(cr2, ry(i+1), n(i+1)*ry(i+2));
            v = reshape(v, r*ry(i+1), d2);
            v = reshape(v.', d2*r, ry(i+1));
            v = v*cr2; # size r+radd, n2, r3
    
            ry(i+1) = r;
    
            u = reshape(u, ry(i), n(i), r);
            v = reshape(v, d2, r, n(i+1), ry(i+2));
    
            # Stuff back
            cry{i} = u;
            cry{i+1} = v;
            
            # Recompute left interface matrices
            # Interface matrix for Y
            Ry{i+1} = Ry{i}*reshape(u, ry(i), n(i)*ry(i+1));
            Ry{i+1} = reshape(Ry{i+1}, ry(i)*n(i), ry(i+1));
            curind = maxvol2(Ry{i+1},'qr',do_qr);
            Ry{i+1} = Ry{i+1}(curind, :);
            # Interface matrices for X
            for j=1:nx
                Rx{i+1,j} = reshape(crX{i,j}, rx(i,j), n(i)*rx(i+1,j));
                Rx{i+1,j} = Rx{i,j}*Rx{i+1,j};
                Rx{i+1,j} = reshape(Rx{i+1,j}, ry(i)*n(i), rx(i+1,j));
                Rx{i+1,j} = Rx{i+1,j}(curind, :);
            end;
        elseif (dir<0)&&(i>1) # right-to-left
            u = u(:,1:r)*diag(s(1:r));
            v = conj(v(:,1:r));
            # kick
            radd = 0; rv = 1;
            if (~last_sweep)&&(kickrank>0)
                if (strcmp(kicktype, 'amr-two'))
                    # Compute the X superblocks
                    ind2 = unique(ceil(rand(ry(i), 1)*(ry(i-1)*n(i-1))));
                    rkick = numel(ind2);
                    curbl = zeros(rkick*n(i)*ry(i+1), nx);
                    for j=1:nx
                        cr1 = reshape(crX{i,j}, rx(i,j)*n(i), rx(i+1,j));
                        cr1 = cr1*Rx{i+1,j};
                        cr1 = reshape(cr1, rx(i,j), n(i)*ry(i+1));
                        cr2 = reshape(crX{i-1,j}, rx(i-1,j), n(i-1)*rx(i,j));                    
                        cr2 = Rx{i-1,j}*cr2; # now its size rx
                        cr2 = reshape(cr2, ry(i-1)*n(i-1), rx(i,j));
                        cr2 = cr2(ind2, :);
                        curbl(:,j) = reshape(cr2*cr1, rkick*n(i)*ry(i+1), 1);
                    end;
                    # Call the function
                    uk = funs(curbl); # rnr x d2
                    uk = reshape(uk, rkick*n(i)*ry(i+1), d2);
                    uk = reshape(uk.', d2*rkick*n(i), ry(i+1));
                    uk = uk / Ry{i+1};
                    uk = reshape(uk, d2*rkick, n(i)*ry(i+1));
                    if (strcmp(pcatype, "svd"))
                        [vk,sk,uk]=svd(uk, "econ");
                        uk = uk(:,1:min(kickrank, size(uk,2)));
                    else
                        uk = uchol(uk, kickrank+1);
                        uk = uk(:,end:-1:max(end-kickrank+1,1));
                    end;                
                else
                    uk = rand(n(i)*ry(i+1), kickrank);
                end;            
    #             uk = rand(n(i)*ry(i+1), kickrank);
                [v,rv]=qr([v,uk], 0);
                radd = size(uk,2);
            end;
            u = [u, zeros(d2*ry(i), radd)];
            u = u*(rv.');
            r = size(v,2);
            
            cr2 = cry{i-1};
            cr2 = reshape(cr2, ry(i-1)*n(i-1), ry(i));
            u = reshape(u, d2, ry(i)*r);
            u = reshape(u.', ry(i), r*d2);
            u = cr2*u;
                   
            u = reshape(u, ry(i-1)*n(i-1)*r, d2);
            u = reshape(u.', d2, ry(i-1), n(i-1), r);
            v = reshape(v.', r, n(i), ry(i+1));
                    
            # Stuff back
            ry(i) = r;
            cry{i-1} = u;
            cry{i} = v;
            
            # Recompute left interface matrices
            # Interface matrix for Y
            Ry{i} = reshape(v, ry(i)*n(i), ry(i+1))*Ry{i+1};
            Ry{i} = reshape(Ry{i}, ry(i), n(i)*ry(i+1));
            curind = maxvol2(Ry{i}.',"qr",do_qr);
            Ry{i} = Ry{i}(:, curind);
            # Interface matrices for X
            for j=1:nx
                Rx{i,j} = reshape(crX{i,j}, rx(i,j)*n(i), rx(i+1,j));
                Rx{i,j} = Rx{i,j}*Rx{i+1,j};
                Rx{i,j} = reshape(Rx{i,j}, rx(i,j), n(i)*ry(i+1));
                Rx{i,j} = Rx{i,j}(:, curind);
            end;
        elseif ((dir>0)&&(i==d))
            # Just stuff back the last core
            newy = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
            newy = reshape(newy, ry(i)*n(i)*ry(i+1), d2);
            cry{i} = reshape(newy.', d2, ry(i), n(i), ry(i+1));
        elseif ((dir<0)&&(i==1))
            # Just stuff back the last core
            newy = u(:,1:r)*diag(s(1:r))*(v(:,1:r)');
            newy = reshape(newy, d2, ry(i), n(i), ry(i+1));
            cry{i} = newy;        
        end;
        
        
        i = i+dir;
        # Reversing, residue check, etc
        cur_order(order_index) = cur_order(order_index) - dir;
        # New direction
        if (cur_order(order_index)==0)
            order_index = order_index+1;
    
            if (verb>0)
                fprintf("=multifuncrs= sweep %d{%d}, max_dy: %3.3e, erank: %g\n", swp, order_index-1, max_dy, sqrt(ry(1:d)'*(n.*ry(2:d+1))/sum(n)));
            end;
    
            if (last_sweep)
                break;
            end;
    
            if (max_dy<eps)&&(dir<0)
                last_sweep=true;
                kickrank=0;
            end;
    
            if (order_index>numel(cur_order)) % New global sweep
                cur_order = block_order;
                order_index = 1;
                #residue
                if (last_sweep)
                    cur_order = d-1;
                end;
    
                max_dy = 0;
                swp = swp+1;
            end;
    
            dir = sign(cur_order(order_index));
            i = i+dir;
        end;
    end
    
    cry{d} = permute(cry{d}, [2,3,1]); # d2 is r(d+1)
    y = cell2core(y, cry);
    
    end
    