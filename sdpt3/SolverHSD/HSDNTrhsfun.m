%%*******************************************************************
%% HSDNTrhsfun: compute the right-hand side vector of the 
%%              Schur complement equation for the NT direction. 
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

    function [rhs,EinvRc,hRd] = HSDNTrhsfun(blk,At,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ);     
    global spdensity 

    m = par.m;   
    if (nargin > 8) 
       corrector = 1; 
    else 
       corrector = 0; 
       hRd = zeros(m+2,1); 
    end       
    hEinvRc = zeros(m+2,1); 
    EinvRc  = cell(size(blk,1),1); 
    if length(sigmu)==1; sigmu = sigmu*ones(1,size(blk,1)); end    
%%
    for p = 1:size(blk,1)
        pblk = blk(p,:); 
        n = sum(pblk{2});  numblk = length(pblk{2});  
        if strcmp(pblk{1},'l')
           if (corrector)
              Rq = dX{p}.*dZ{p}; 
           else
              Rq = sparse(n,1); 
              tmp  = par.dd{p}.*Rd{p};
              tmp2 = mexMatvec(At{p},tmp,1);
              hRd = hRd + tmp2;
           end
           EinvRc{p} = sigmu(p)./Z{p}-X{p} -Rq./Z{p};
           tmp2 = mexMatvec(At{p},EinvRc{p},1);  
           hEinvRc = hEinvRc + tmp2;
	elseif strcmp(pblk{1},'q') 
           w = sqrt(par.gamz{p}./par.gamx{p}); 
   	   if (corrector)
              hdx = qops(pblk,w,par.ff{p},5,dX{p}); 
              hdz = qops(pblk,w,par.ff{p},6,dZ{p}); 
              hdxdz = Arrow(pblk,hdx,hdz);
              vv = qops(pblk,w,par.ff{p},5,X{p}); 
              Vihdxdz = Arrow(pblk,vv,hdxdz,1); 
              Rq = qops(pblk,w,par.ff{p},6,Vihdxdz); 
           else
              Rq = sparse(n,1); 
              tmp  = par.dd{p}.*Rd{p} + qops(pblk,qops(pblk,Rd{p},par.ee{p},1),par.ee{p},3);
              tmp2 = mexMatvec(At{p},tmp,1);
              hRd = hRd + tmp2;
           end
           EinvRc{p} = qops(pblk,-sigmu(p)./(par.gamz{p}.*par.gamz{p}),Z{p},4)-X{p}-Rq;
           tmp2 = mexMatvec(At{p},EinvRc{p},1);         
           hEinvRc = hEinvRc + tmp2;
        elseif strcmp(pblk{1},'s') 
           n2 = pblk{2}.*(pblk{2}+1)/2; 
           if (corrector)
              hdZ = Prod3(pblk,par.G{p},dZ{p},par.G{p}',1); 
              hdX = spdiags(-par.sv{p},0,n,n)-hdZ;          
              tmp = Prod2(pblk,hdX,hdZ,0);  
              tmp = 0.5*(tmp+tmp');
              if (numblk == 1) 
                 d = par.sv{p};
       	         e = ones(pblk{2},1); 
     	         Rq = 2*tmp./(d*e'+e*d'); 
                 if (nnz(Rq) <= spdensity*n2); Rq = sparse(Rq); end
              else
                 Rq = sparse(n,n);
                 s = [0, cumsum(pblk{2})]; 
                 for i = 1:numblk
                     pos = [s(i)+1 : s(i+1)]; 
                     d = par.sv{p}(pos); e = ones(length(pos),1); 
                     Rq(pos,pos) = 2*tmp(pos,pos)./(d*e' + e*d'); 
                 end
              end
           else
              Rq = sparse(n,n); 
              EinvRc{p} = Prod3(pblk,par.W{p},Rd{p},par.W{p},1,par.nzlistAy{p});
              tmp2 = AXfun(pblk,At(p,:),par.permA(p,:),EinvRc(p)); 
              hRd = hRd + tmp2;
           end 
           tmp = spdiags(sigmu(p)./par.sv{p} -par.sv{p},0,n,n);
           EinvRc{p} = Prod3(pblk,par.G{p}',tmp-Rq,par.G{p},1);
           tmp2 = AXfun(pblk,At(p,:),par.permA(p,:),EinvRc(p)); 
           hEinvRc = hEinvRc + tmp2;  
        end 
    end
%% 
    rhs = rp + hRd - hEinvRc; 
    rhs(m+1) = rhs(m+1) + (par.mu/par.tau - par.kap); 
    if (corrector)
       rhs(m+1) = rhs(m+1) - par.dtau*par.dkap/par.tau; 
    end
%%*******************************************************************
