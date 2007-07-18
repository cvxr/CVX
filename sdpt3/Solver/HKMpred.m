%%*******************************************************************
%% HKMpred: Compute (dX,dy,dZ) for the H..K..M direction. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

   function [par,dX,dy,dZ,coeff,L,hRd] = ...
           HKMpred(blk,At,par,rp,Rd,sigmu,X,Z,invZchol);
       
    global schurfun schurfun_par
%%
%% compute HKM scaling 
%%
   Zinv = cell(size(blk,1),1);   dd = cell(size(blk,1),1);
   gamx = cell(size(blk,1),1);   gamz = cell(size(blk,1),1); 
   for p = 1:size(blk,1) 
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      numblk = length(pblk{2}); 
      if strcmp(pblk{1},'l')
         Zinv{p} = 1./Z{p};
         dd{p} = X{p}./Z{p}; %% do not add perturbation, it badly affects cre-a
      elseif strcmp(pblk{1},'q')
         gaptmp  = qops(pblk,X{p},Z{p},1); 
         gamz2   = qops(pblk,Z{p},Z{p},2); 
         gamz{p} = sqrt(gamz2); 
         Zinv{p} = qops(pblk,-1./gamz2,Z{p},4);  
         dd{p}   = qops(pblk,gaptmp./gamz2,ones(n,1),4); 
      elseif strcmp(pblk{1},'s')
         if (numblk == 1)
            Zinv{p} = Prod2(pblk,full(invZchol{p}),invZchol{p}',1);
            %% to fix the anonmaly when Zinv has very small elements 
            if (par.iter==2 | par.iter==3) & ~issparse(Zinv{p}); 
               Zinv{p} = Zinv{p} + 1e-16; 
            end
         else
            Zinv{p} = Prod2(pblk,invZchol{p},invZchol{p}',1);            
         end
      end
   end
   par.Zinv = Zinv; par.gamx = gamx; par.gamz = gamz; par.dd = dd; 
%%
%% compute schur matrix
%%
    m = length(rp);
    schur = sparse(m,m); 
    UU = [];  EE = [];  Afree = [];  
    dX = cell(size(blk,1),1); dy = []; dZ = cell(size(blk,1),1); 
%%
    for p = 1:size(blk,1)
       pblk = blk(p,:); 
       if strcmp(pblk{1},'l')
          [schur,UU,EE] = schurmat_lblk(blk,At,par,schur,UU,EE,p,par.dd);
       elseif strcmp(pblk{1},'q');  
          [schur,UU,EE] = schurmat_qblk(blk,At,par,schur,UU,EE,p,par.dd,par.Zinv,X);
       elseif strcmp(pblk{1},'s')
          if isempty(schurfun{p})
             schur = schurmat_sblk(blk,At,par,schur,p,X,par.Zinv); 
          elseif isstr(schurfun{p}) 
             schurtmp = sparse(m,m);
             if ~isempty(par.permZ{p})
                Zpinv = Zinv{p}(par.permZ{p},par.permZ{p}); 
                Xp = X{p}(par.permZ{p},par.permZ{p}); 
             else
                Xp = X{p}; 
                Zpinv = Zinv{p};
             end
             eval(['schurtmp = ',schurfun{p},'(Xp,Zpinv,schurfun_par(p,:));']); 
             schur = schur + schurtmp;
          end
       elseif strcmp(pblk{1},'u') 
          Afree = [Afree, At{p}'];
       end
    end
%%
%% compute rhs 
%%
    [rhs,EinvRc,hRd] = HKMrhsfun(blk,At,par,X,Z,rp,Rd,sigmu);
%%
%% solve linear system
%%    
    [xx,coeff,L] = linsysolve(par,schur,UU,Afree,EE,rhs); 
%% 
%% compute (dX,dZ)
%%  
    [dX,dy,dZ] = HKMdirfun(blk,At,par,Rd,EinvRc,X,xx,m);
%%*******************************************************************
