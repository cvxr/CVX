%%**********************************************************************
%% NTpred: Compute (dX,dy,dZ) for NT direction. 
%%                       
%% compute SVD of Xchol*Zchol via eigenvalue decompostion of
%%     Zchol * X * Zchol' = V * diag(sv2) * V'. 
%% compute W satisfying W*Z*W = X. 
%%     W = G'*G,  where G = diag(sqrt(sv)) * (invZchol*V)'
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**********************************************************************

 function [par,dX,dy,dZ,coeff,L,hRd] = ...
          NTpred(blk,At,par,rp,Rd,sigmu,X,Z,Zchol,invZchol);

    global schurfun schurfun_par 
%%
%% compute NT scaling matrix
%%
    [par.W,par.G,par.sv,par.gamx,par.gamz,par.dd,par.ee,par.ff] = ...
     NTscaling(blk,X,Z,Zchol,invZchol);
%%
%% compute schur matrix
%%
    m = length(rp); 
    schur = sparse(m,m); 
    UU = []; EE = []; Afree = []; 
    dX = cell(size(blk,1),1); dy = []; dZ = cell(size(blk,1),1); 
%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       if strcmp(pblk{1},'l')
          [schur,UU,EE] = schurmat_lblk(blk,At,par,schur,UU,EE,p,par.dd);
       elseif strcmp(pblk{1},'q');       
          [schur,UU,EE] = schurmat_qblk(blk,At,par,schur,UU,EE,p,par.dd,par.ee);
       elseif strcmp(pblk{1},'s')
          if isempty(schurfun{p})
             schur = schurmat_sblk(blk,At,par,schur,p,par.W); 
          elseif isstr(schurfun{p}) 
             schurtmp = sparse(m,m);
             if ~isempty(par.permZ{p})
                Wp = par.W{p}(par.permZ{p},par.permZ{p}); 
             else
                Wp = par.W{p};
             end
             eval(['schurtmp = ',schurfun{p},'(Wp,Wp,schurfun_par(p,:));']); 
             schur = schur + schurtmp;
          end
       elseif strcmp(pblk{1},'u')            
          Afree = [Afree, At{p}'];
       end
    end
%%
%% compute rhs
%%
    [rhs,EinvRc,hRd] = NTrhsfun(blk,At,par,X,Z,rp,Rd,sigmu);
%%
%% solve linear system
%%
    [xx,coeff,L] = linsysolve(par,schur,UU,Afree,EE,rhs); 
%%
%% compute (dX,dZ)
%%
    [dX,dy,dZ] = NTdirfun(blk,At,par,Rd,EinvRc,xx,m); 
%%**********************************************************************
