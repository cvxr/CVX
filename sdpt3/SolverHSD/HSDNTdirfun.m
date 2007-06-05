%%*******************************************************************
%% HSDNTdirfun: compute (dX,dZ), given dy, for the NT direction.
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

    function [par,dX,dy,dZ] = HSDNTdirfun(blk,At,par,Rd,EinvRc,xx); 
    
    global solve_ok

    dX = cell(size(blk,1),1); dZ = cell(size(blk,1),1); dy = [];
    if (any(isnan(xx)) | any(isinf(xx)))
       solve_ok = 0;
       fprintf('\n  HSDNTdirfun: solution contains NaN or inf.');
       return;
    end
%%
    m = par.m; 
    dy2 = xx(1:m+2); 
%%
    for p=1:size(blk,1)
       pblk = blk(p,:);  
       if strcmp(pblk{1},'l')
          dZ(p) = ops(Rd(p),'-',Atyfun(pblk,At(p,:),[],[],dy2));
          tmp   = par.dd{p}.*dZ{p};
          dX{p} = EinvRc{p} - tmp;
       elseif strcmp(pblk{1},'q')
          dZ(p) = ops(Rd(p),'-',Atyfun(pblk,At(p,:),[],[],dy2));
          tmp = par.dd{p}.*dZ{p} + qops(pblk,qops(pblk,dZ{p},par.ee{p},1),par.ee{p},3); 
          dX{p} = EinvRc{p} - tmp;       
       elseif strcmp(pblk{1},'s') 
          dZ(p) = ops(Rd(p),'-',Atyfun(pblk,At(p,:),par.permA(p,:),par.isspAy(p),dy2)); 
          tmp   = Prod3(pblk,par.W{p},dZ{p},par.W{p},1); 
          dX{p} = EinvRc{p}-tmp;
       end
    end 
    dy  = dy2(1:m); 
    par.dtau = dy2(m+1); 
    par.dtheta = dy2(m+2);
    par.dkap = (par.mu./par.tau - par.kap) - par.kap*(par.dtau/par.tau); 
%%*******************************************************************
