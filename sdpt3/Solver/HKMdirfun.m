%%*******************************************************************
%% HKMdirfun: compute (dX,dZ), given dy, for the HKM direction.
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

    function [dX,dy,dZ] = HKMdirfun(blk,At,par,Rd,EinvRc,X,xx,m); 

    global solve_ok

    dX = cell(size(blk,1),1); dZ = cell(size(blk,1),1); dy = [];
    if (any(isnan(xx)) | any(isinf(xx)))
       solve_ok = 0;
       fprintf('\n  linsysolve: solution contains NaN or inf');
       return;
    end
%%
    dy = xx(1:m); 
    count = m;
%%
    for p=1:size(blk,1)
       pblk = blk(p,:);  
       if strcmp(pblk{1},'l')
          %%dZ{p} = Rd{p} - At{p}*dy; 
          dZ(p) = ops(Rd(p),'-',Atyfun(pblk,At(p,:),[],[],dy)); 
          dX{p} = EinvRc{p} - par.dd{p}.*dZ{p};
       elseif strcmp(pblk{1},'q') 
          %%dZ{p} = Rd{p} - At{p}*dy;  
          dZ(p) = ops(Rd(p),'-',Atyfun(pblk,At(p,:),[],[],dy)); 
          tmp = par.dd{p}.*dZ{p} ...
                + qops(pblk,qops(pblk,dZ{p},par.Zinv{p},1),X{p},3) ...
                + qops(pblk,qops(pblk,dZ{p},X{p},1),par.Zinv{p},3); 
          dX{p} = EinvRc{p} - tmp;
       elseif strcmp(pblk{1},'s')
          %%dZ{p} = Rd{p} -smat(pblk,At{p}*dy(par.permA(p,:)),par.isspAy(p));
          dZ(p) = ops(Rd(p),'-',Atyfun(pblk,At(p,:),par.permA(p,:),par.isspAy(p),dy));
          tmp = Prod3(pblk,X{p},dZ{p},par.Zinv{p},0); 
          tmp = 0.5*(tmp+tmp'); 
          dX{p} = EinvRc{p}-tmp;
       elseif strcmp(pblk{1},'u'); 
          n = sum(pblk{2}); 
          dZ{p} = zeros(n,1); 
          dX{p} = xx(count+[1:n]); 
          count = count + n;
       end
    end
%%*******************************************************************
