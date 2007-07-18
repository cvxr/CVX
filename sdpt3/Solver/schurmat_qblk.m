%%*******************************************************************
%% schurmat_qblk: compute schur matrix corresponding to SOCP blocks.
%%
%% HKM direction: output = schur + Ax*Ae' + Ae*Ax' - Ad*Ad'
%% NT  direction: output = schur + Ae*Ae' - Ad*Ad'
%%
%% where schur = A*D*A', and Ad is the modification to ADA' 
%% so that the latter is positive definite. 
%%
%% [schur,UU,EE] = schurmat_qblk(blk,At,schur,UU,EE,p,dd,ee,xx);
%% 
%% UU: stores the dense columns of Ax, Ae, Ad, and possibly 
%%     those of A*D^{1/2}. It has the form UU = [Ax Ae Ad]. 
%% EE: stores the assocaited (2,2) block matrix when the
%%     output matrix is expressed as an augmented matrix.
%%     It has the form EE = [0 -lam 0; -lam 0 0; 0 0 I].
%%
%% options = 0, HKM
%%         = 1, NT
%% 
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

   function [schur,UU,EE] = schurmat_qblk(blk,At,par,schur,UU,EE,p,dd,ee,xx);
   
   global idxdenAq nnzschur_qblk

   if (nargin == 10); options = 0; else; options = 1; end; 
   iter = par.iter; 
      
   if isempty(EE) 
      count = 0; 
   else 
      count = max(max(EE(:,2)),max(EE(:,1))); 
   end
   pblk = blk(p,:); n = sum(pblk{2}); numblk = length(pblk{2}); 
%%   
   Ae = qprod(pblk,At{p}',ee{p}); 
   if (options == 0) 
      Ax = qprod(pblk,At{p}',xx{p}); 
   end
   idxden = checkdense(Ae);
   ddsch = dd{p};    
   if ~isempty(idxden); 
      spcolidx = setdiff([1:numblk],idxden); 
      s = 1 + [0, cumsum(pblk{2})];
      idx = s(idxden); 
      tmp = zeros(n,1); 
      tmp(idx) = sqrt(2*abs(ddsch(idx))); 
      Ad = qprod(pblk,At{p}',tmp); 
      ddsch(idx) = abs(ddsch(idx)); 
      if (options == 0) 
         len = length(idxden); 
         gamzsub = par.gamz{p}(idxden);
         lam = gamzsub.*gamzsub;
         UU = [UU, Ax(:,idxden), Ae(:,idxden)*spdiags(lam,0,len,len), Ad(:,idxden)]; 
         tmp = count+[1:len]'; 
         EE = [EE; [tmp, len+tmp, -lam; len+tmp, tmp, -lam; ...
                    2*len+tmp, 2*len+tmp, ones(len,1)] ];
         count = count+3*len;  
         Ax = Ax(:,spcolidx); Ae = Ae(:,spcolidx); 
         tmp = Ax*Ae'; 
         schur = schur + (tmp + tmp');
      else
         len = length(idxden);
         w2 = par.gamz{p}./par.gamx{p}; 
         lam = w2(idxden); 
         UU = [UU, Ae(:,idxden)*spdiags(sqrt(lam),0,len,len), Ad(:,idxden)]; 
         tmp = count+[1:len]'; 
         EE = [EE; [tmp, tmp, -lam; len+tmp, len+tmp, ones(len,1)] ]; 
         count = count + 2*len; 
         Ae = Ae(:,spcolidx);      
         schur = schur + Ae*Ae';
      end
   else
      if (options == 0)
         tmp = Ax*Ae'; 
         schur = schur + (tmp+tmp');
      else 
         tmp = Ae*Ae'; 
         schur = schur + tmp; 
      end
   end
   if (iter==1)
      idxdenAq{p} = checkdense(At{p}'); 
   end
   if ~isempty(idxdenAq{p});
      idxden = idxdenAq{p};  
      len = length(idxden);               
      Ad = At{p}(idxden,:)'*spdiags(sqrt(abs(ddsch(idxden))),0,len,len); 
      UU = [UU, Ad];
      tmp = count+[1:len]'; 
      EE = [EE; [tmp, tmp, -sign(ddsch(idxden))]]; 
      count = count + len; 
      ddsch(idxden) = zeros(len,1); 
   end  
   schurtmp = At{p}' *spdiags(ddsch,0,n,n) *At{p}; 
   schur = schur + schurtmp;
%%*******************************************************************
