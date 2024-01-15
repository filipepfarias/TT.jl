
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
            ps=tt.ps;
            core=tt.core;
            n=tt.n;
            r=tt.r;
            
            ry=y.r;
            psy=y.ps;
            cry=y.core;
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
                cr=reshape(cr,ry[i+1],n[i+1]*ry[i+2]);
                cr=transpose(cr);
                cr,rm=qr(cr); cr = Array(cr); 
                ry[i+1]=size(cr,2);
                ind=maxvol2(cr); 
                r1=cr[ind,:];
                cr=cr/r1; #cr[abs.(cr) .< eps() ] .= 0.0  ################################ Needs to be checked.
                cr=transpose(cr); 
                cry[pos1:pos1+ry[i+1]*n[i+1]*ry[i+2]-1]=cr[:];
                pos1=pos1-ry[i]*n[i]*ry[i+1];
                cr2=cr2*transpose(r1*rm); 
                cry[pos1:pos1+ry[i]*n[i]*ry[i+1]-1]=cr2[:];
                #Take phi matrix; convolve from right with current cors of V
                cr0=core[ps[i+1]:ps[i+2]-1];
                cr0=reshape(cr0,r[i+1]*n[i+1],r[i+2]); 
                cr0=cr0*phi[i+2]; #cr0 is now r(i)*n(i)*ry(i+1);
                cr0=reshape(cr0,r[i+1],n[i+1]*ry[i+2]);
                phi[i+1]=cr0[:,ind]; 
                #  psy=cumsum([1;n.*ry(1:d).*ry(2:d+1)]);
                #tmpy_core=cry;
                #tmpy_r=ry;
                #tmpy_ps=psy;
                #keyboard
                
                #Orthogonalization & maxvol is "local operation" (touches two
                #cores)
                #Computation of elements requries storage of phi matrices   
            end
            
            #Truncate cry
            #pos1=pos1-n(1)*r(2);
            cry=cry[pos1:prod(size(cry))];
            tmpy_core=cry;
            tmpy_r=ry;
            tmpy_ps=psy;
            #keyboard
            swp=1;
            
            yold=[];
            not_converged = true;
            pos1=1;
            last_sweep = false;
            while ( swp < n_sweeps && not_converged )
                max_er=0;
                cry_old=cry;
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
                    display(i);display(core);
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
                    s=s(1:r2); u=u(:,1:r2); v=v(:,1:r2);
                    
                    v=conj(v)*diag(s);
                    
                    if ((swp % dropsweeps)!=0 && (swp > 1) && !last_sweep)
                        u = [cry1 u];
                        v = [transpose(cry2) v];
                        u,rv=qr(u);
                        v = v*transpose(rv);         
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
                    u=u/r1; v=v*r1'; v=transpose(v);
                    #Compute new phi
                    ry[i+1]=r2;
                    cry[pos1:pos1+ry[i]*n[i]*ry[i+1]-1]=u[:];
                    pos1=pos1+ry[i]*n[i]*ry[i+1];
                    cry[pos1:pos1+ry[i+1]*n[i+1]*ry[i+2]-1]=v[:];
                    
                    #Now we have to: cr1 with phi from the left 
                    phi[i+1]=cr1[ind,:]; #phi[i+1]=phi[i+1];    
                end
                
                #Truncate local memory
                cry=cry[1:pos1+ry[d]*n[d]*ry[d+1]-1];
                psy=cumsum([1;n.*ry[1:d].*ry[2:d+1]]);
                tmpy_core=cry;
                tmpy_r=ry;
                tmpy_ps=psy;
                #keyboard;
                cry_old=cry;
                #m=numel(cry);
                #cry=[zeros(size(cry)),cry];
                #pos1=pos1+m;
                #cry2=cry_old(psy(i+1):psy(i+2)-1);
                cry=cry[psy[d]:psy[d+1]-1]; #Start--only two cores
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
                    cr=fun(cr); #Elements are evaluated here!
                    cr=reshape(cr,ry[i]*n[i],n[i+1]*ry[i+2]);
                    #Check for local approximation of cr for the error
                    #cry1=cry(pos1-ry(i)*n(i)*ry(i+1):pos1-1);
                    cry1=cry_old[psy[i]:psy[i+1]-1];
                    #cry2=cry(ry(:ry(i+1)*n(i+1)*ry(i+2));
                    #cry1=cry(1:ry(i)*n(i)*ry(i+1));
                    #cry2=cry(ry(i)*n(i)*ry(i+1)+1:ry(i)*n(i)*ry(i+1)+ry(i+1)*n(i+1)*ry(i+2));
                    cry2=cry[1:ry[i+1]*n[i+1]*ry[i+2]];
                    cry[1:ry[i+1]*n[i+1]*ry[i+2]]=[];
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
                        v = [transpose(cry2) v];
                        v,rv=qr(v);
                        u = u*transpose(rv);
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
                    v=v/r1; v=transpose(v);
                    u=u*r1';
                    #Compute new phi;
                    ry[i+1]=r2;
                    u=u[:]; u=u'; v=v[:]; v=v';
                    cry=[u v cry];
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
                    
                end
                #Truncate local memory
                #cry=cry(pos1:numel(cry));
                
                # keyboard;
                psy=cumsum([1;n.*ry[1:d].*ry[2:d+1]]);
                
                tmpy_core=cry;
                tmpy_r=ry;
                tmpy_ps=psy;
                #keyboard;
                if ( isempty(yold) )
                    yold=y;
                    er_nrm=1;
                else
                    
                    er_nrm=norm(yold-y)/norm(y);
                    yold=y;
                end
                if ( verb )
                    @sprintf("sweep=%d, er=%3.2e er_nrm=%3.2e \n",swp,max_er,er_nrm);
                end
                if (last_sweep)
                    break;
                end;
                if (er_nrm<precision)
                    last_sweep=true;
                end;
                swp=swp+1;
            end     
            psy=cumsum([1;n.*ry[1:d].*ry[2:d+1]]);
            
            tmpy_core=cry;
            tmpy_r=ry;
            tmpy_ps=psy;
            
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
                return ind
            end
            
            if ( do_qr ) 
                a,_ = qr(a);
            end
            
            
            #Initialize
            if ( do_lu_full )
                ind = lu_full(a); 
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
                k = ind[big_ind[2]]; #This is the current row that we are using
                b = b + b[:,big_ind[2]]*(b[k,:]-b[big_ind[1],:])'/b[big_ind];
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
                reort_flag=~isempty(find(norm_unew <= 0.25*norm_uadd,1));
                
                unew,_=qr(unew); #Here it is ok.
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
                disp("Reort failed to reort!");
                y,_ =qr([u unew]);
            else
                y=[u unew];
            end
            #if ( ~isempty(z) ) 
            #   y = y - 2*z*(z'*y);
            #end
            return y
        end