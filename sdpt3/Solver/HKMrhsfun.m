%%*******************************************************************
%% HKMrhsfun: compute the right-hand side vector of the 
%%            Schur complement equation for the  HKM direction. 
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

  function [rhs,EinvRc,hRd] = HKMrhsfun(blk,At,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ)

    m = length(rp);   
    if (nargin > 8)
       corrector = 1;  
    else
       corrector = 0; 
       hRd = zeros(m,1);
    end       
    hEinvRc = zeros(m,1); 
    EinvRc = cell(size(blk,1),1); 
    rhsfree = []; 
%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       n = sum(pblk{2}); 
       if strcmp(pblk{1},'l')
	  if iscell(sigmu)
             EinvRc{p} = sigmu{p}./Z{p} -X{p};
	  else
             EinvRc{p} = sigmu./Z{p} -X{p};
          end
          Rq = sparse(n,1); 
          if (corrector) & (norm(par.parbarrier{p})==0)
             Rq = dX{p}.*dZ{p}./Z{p};
          else
             tmp  = par.dd{p}.*Rd{p};
             tmp2 = mexMatvec(At{p},tmp,1);
             hRd = hRd + tmp2;
          end
	  EinvRc{p} = EinvRc{p} - Rq; 
          tmp2 = mexMatvec(At{p,1},EinvRc{p},1);  
          hEinvRc = hEinvRc + tmp2;
       elseif strcmp(pblk{1},'q') 
	  if iscell(sigmu)
             EinvRc{p} = qops(pblk,sigmu{p},par.Zinv{p},3) -X{p};
	  else
             EinvRc{p} = sigmu*par.Zinv{p} -X{p};
          end
          Rq = sparse(n,1); 
          if (corrector) & (norm(par.parbarrier{p})==0)
  	     ff{p} = qops(pblk,1./par.gamz{p},Z{p},3);
             hdx = qops(pblk,par.gamz{p},ff{p},5,dX{p}); 
             hdz = qops(pblk,par.gamz{p},ff{p},6,dZ{p}); 
             hdxdz = Arrow(pblk,hdx,hdz);
             Rq = qops(pblk,par.gamz{p},ff{p},6,hdxdz); 
          else
             tmp = par.dd{p}.*Rd{p} ...
                   + qops(pblk,qops(pblk,Rd{p},par.Zinv{p},1),X{p},3) ...
                   + qops(pblk,qops(pblk,Rd{p},X{p},1),par.Zinv{p},3); 
             tmp2 = mexMatvec(At{p,1},tmp,1);
             hRd = hRd + tmp2;
          end
	  EinvRc{p} = EinvRc{p} - Rq; 
          tmp2 = mexMatvec(At{p,1},EinvRc{p},1);         
          hEinvRc = hEinvRc + tmp2;
       elseif strcmp(pblk{1},'s') 
	  if iscell(sigmu)	    
             %%ss = [0,cumsum(pblk{2})]; 
             %%sigmuvec = zeros(n,1); 
             %%for k = 1:length(pblk{2}); 
             %%   sigmuvec(ss(k)+1:ss(k+1)) = sigmu{p}(k)*ones(pblk{2}(k),1); 
             %%end
             sigmuvec = mexexpand(pblk{2},sigmu{p}); 
             EinvRc{p} = par.Zinv{p}*spdiags(sigmuvec,0,n,n) -X{p};
          else
             EinvRc{p} = sigmu*par.Zinv{p} -X{p};
          end
          Rq = sparse(n,n);
          if (corrector) & (norm(par.parbarrier{p})==0)
             Rq = Prod3(pblk,dX{p},dZ{p},par.Zinv{p},0); 
	     Rq = 0.5*(Rq+Rq');
          else
             tmp = Prod3(pblk,X{p},Rd{p},par.Zinv{p},0,par.nzlistAy{p}); 
             tmp = 0.5*(tmp+tmp'); 
             tmp2 = AXfun(pblk,At(p,:),par.permA(p,:),{tmp}); 
             hRd = hRd + tmp2;
          end 
	  EinvRc{p} = EinvRc{p} - Rq; 
          tmp2 = AXfun(pblk,At(p,:),par.permA(p,:),EinvRc(p)); 
          hEinvRc = hEinvRc + tmp2;
       elseif strcmp(pblk{1},'u') 
          rhsfree = [rhsfree; Rd{p}]; 
       end
    end
%%
    rhs = rp + hRd - hEinvRc; 
    rhs = full([rhs; rhsfree]);  
%%*******************************************************************
