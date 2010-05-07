%%*****************************************************************
%% HKMcorr: corrector step for the HKM direction. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*****************************************************************

  function [dX,dy,dZ,resnrm,EinvRc] = HKMcorr(blk,At,par,rp,Rd,sigmu,hRd,...
            dX,dZ,coeff,L,X,Z);

    global matfct_options solve_ok 

    printlevel = par.printlevel;
%%
    [rhs,EinvRc] = HKMrhsfun(blk,At,par,X,Z,rp,Rd,sigmu,hRd,dX,dZ);
    m = length(rp); ncolU = size(coeff.mat12,2); 
    rhs = [rhs; zeros(m+ncolU-length(rhs),1)];
%%
    solve_ok = 1; resnrm = norm(rhs);
    if strcmp(matfct_options,'chol') | strcmp(matfct_options,'spchol') ...
       | strcmp(matfct_options,'ldl') | strcmp(matfct_options,'spldl')
       [xx,resnrm,solve_ok] = symqmr(coeff,rhs,L,[],[],printlevel);
       if (solve_ok<=0) & (printlevel)
          fprintf('\n  warning: symqmr fails: %3.1f.',solve_ok); 
       end
    else
       [xx,resnrm,solve_ok] = mybicgstab(coeff,rhs,L,[],[],printlevel);
       if (solve_ok<=0) & (printlevel)
          fprintf('\n  warning: bicgstab fails: %3.1f.',solve_ok); 
       end
    end
    if (printlevel>=3); fprintf('%2.0d ',length(resnrm)-1); end
%%
    [dX,dy,dZ] = HKMdirfun(blk,At,par,Rd,EinvRc,X,xx,m);
%%*****************************************************************
