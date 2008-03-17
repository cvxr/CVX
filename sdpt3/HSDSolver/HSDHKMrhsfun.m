%%*******************************************************************
%% HSDHKMrhsfun: compute the right-hand side vector of the 
%%               Schur complement equation for the  HKM direction. 
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

  function [rhs,EinvRc,hRd] = HSDHKMrhsfun(blk,At,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ)

    m = par.m;   
    if (nargin > 8)
       corrector = 1; 
    else 
       corrector = 0; 
       hRd = zeros(m+2,1); 
    end       
    hEinvRc = zeros(m+2,1); 
    EinvRc = cell(size(blk,1),1); 
    if length(sigmu)==1; sigmu = sigmu*ones(1,size(blk,1)); end    
%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       n = sum(pblk{2}); 
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
          tmp2 = mexMatvec(At{p,1},EinvRc{p},1);  
          hEinvRc = hEinvRc + tmp2;
       elseif strcmp(pblk{1},'q') 
          if (corrector)
  	     ff{p} = qops(pblk,1./par.gamz{p},Z{p},3);
             hdx = qops(pblk,par.gamz{p},ff{p},5,dX{p}); 
             hdz = qops(pblk,par.gamz{p},ff{p},6,dZ{p}); 
             hdxdz = Arrow(pblk,hdx,hdz);
             Rq = qops(pblk,par.gamz{p},ff{p},6,hdxdz); 
          else
             Rq  = sparse(n,1); 
             tmp = par.dd{p}.*Rd{p} ...
                   + qops(pblk,qops(pblk,Rd{p},par.Zinv{p},1),X{p},3) ...
                   + qops(pblk,qops(pblk,Rd{p},X{p},1),par.Zinv{p},3); 
             tmp2 = mexMatvec(At{p,1},tmp,1);
             hRd = hRd + tmp2;
          end
          EinvRc{p} = sigmu(p)*par.Zinv{p}-X{p} -Rq;
          tmp2 = mexMatvec(At{p,1},EinvRc{p},1);         
          hEinvRc = hEinvRc + tmp2;
       elseif strcmp(pblk{1},'s') 
          if (corrector)
             Rq = Prod3(pblk,dX{p},dZ{p},par.Zinv{p},0); 
	     Rq = 0.5*(Rq+Rq'); 
          else
             Rq = sparse(n,n);
             tmp = Prod3(pblk,X{p},Rd{p},par.Zinv{p},0,par.nzlistAy{p}); 
             EinvRc{p} = 0.5*(tmp+tmp'); 
             tmp2 = AXfun(pblk,At(p,:),par.permA(p,:),EinvRc(p)); 
             hRd = hRd + tmp2;
          end 
          EinvRc{p} = sigmu(p)*par.Zinv{p}-X{p}- Rq;
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
