%%*******************************************************************
%% NTrhsfun: compute the right-hand side vector of the 
%%           Schur complement equation for the NT direction. 
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

    function [rhs,EinvRc,hRd] = NTrhsfun(blk,At,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ);

    spdensity = par.spdensity; 
    m = length(rp);   
    if (nargin > 8) 
       corrector = 1; 
    else 
       corrector = 0; 
       hRd = zeros(m,1); 
    end       
    hEinvRc = zeros(m,1); 
    EinvRc  = cell(size(blk,1),1); 
    rhsfree = []; 
%%
    for p = 1:size(blk,1)
        pblk = blk(p,:); 
        n = sum(pblk{2});  numblk = length(pblk{2});  
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
           tmp2 = mexMatvec(At{p},EinvRc{p},1);  
           hEinvRc = hEinvRc + tmp2;
	elseif strcmp(pblk{1},'q') 
	   if iscell(sigmu)
              EinvRc{p} = qops(pblk,-sigmu{p}./(par.gamz{p}.*par.gamz{p}),Z{p},4) -X{p};
	   else
              EinvRc{p} = qops(pblk,-sigmu./(par.gamz{p}.*par.gamz{p}),Z{p},4) -X{p};
           end
           Rq = sparse(n,1); 
   	   if (corrector) & (norm(par.parbarrier{p})==0)
              w = sqrt(par.gamz{p}./par.gamx{p}); 
              hdx = qops(pblk,w,par.ff{p},5,dX{p}); 
              hdz = qops(pblk,w,par.ff{p},6,dZ{p}); 
              hdxdz = Arrow(pblk,hdx,hdz);
              vv = qops(pblk,w,par.ff{p},5,X{p}); 
              Vihdxdz = Arrow(pblk,vv,hdxdz,1); 
              Rq = qops(pblk,w,par.ff{p},6,Vihdxdz); 
           else
              tmp  = par.dd{p}.*Rd{p} + qops(pblk,qops(pblk,Rd{p},par.ee{p},1),par.ee{p},3);
              tmp2 = mexMatvec(At{p},tmp,1);
              hRd = hRd + tmp2;
           end
	   EinvRc{p} = EinvRc{p} - Rq; 
           tmp2 = mexMatvec(At{p},EinvRc{p},1);         
           hEinvRc = hEinvRc + tmp2;
        elseif strcmp(pblk{1},'s') 
           n2 = pblk{2}.*(pblk{2}+1)/2; 
	   if iscell(sigmu)
     	      %%ss = [0,cumsum(pblk{2})]; 
              %%sigmuvec = zeros(n,1); 
              %%for k = 1:length(pblk{2}); 
              %%   sigmuvec(ss(k)+1:ss(k+1)) = sigmu{p}(k)*ones(pblk{2}(k),1); 
              %%end
              sigmuvec = mexexpand(pblk{2},sigmu{p}); 
              tmp = spdiags(sigmuvec./par.sv{p} -par.sv{p},0,n,n);
           else
              tmp = spdiags(sigmu./par.sv{p} -par.sv{p},0,n,n);
           end
           EinvRc{p} = Prod3(pblk,par.G{p}',tmp,par.G{p},1);
           Rq = sparse(n,n); 
           if (corrector) & (norm(par.parbarrier{p})==0)
              hdZ = Prod3(pblk,par.G{p},dZ{p},par.G{p}',1); 
              hdX = spdiags(qops(pblk,par.parbarrier{p}',1./par.sv{p},3)-par.sv{p},0,n,n)-hdZ; 
              tmp = Prod2(pblk,hdX,hdZ,0);  
              tmp = 0.5*(tmp+tmp');
              if (numblk == 1) 
                 d = par.sv{p};
                 e = ones(pblk{2},1); 
                 Rq = 2*tmp./(d*e'+e*d'); 
                 if (nnz(Rq) <= spdensity*n2); Rq = sparse(Rq); end 
              else
                 Rq = sparse(n,n);
                 ss = [0, cumsum(pblk{2})]; 
                 for i = 1:numblk
                    pos = [ss(i)+1 : ss(i+1)]; 
                    d = par.sv{p}(pos); e = ones(length(pos),1); 
                    Rq(pos,pos) = 2*tmp(pos,pos)./(d*e' + e*d'); 
                 end
              end
              Rq = Prod3(pblk,par.G{p}',Rq,par.G{p},1);
           else
              tmp = Prod3(pblk,par.W{p},Rd{p},par.W{p},1,par.nzlistAy{p});
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
