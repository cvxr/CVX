%%*******************************************************************
%% detect_lblk: detect diagonal blocks in the SDP data. 
%%
%% [blk,At,C,diagblkinfo,blockchange,parbarrier,X,Z] = ...
%%           detect_lblk(blk,At,C,b,parbarrier,X,Z); 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%******************************************************************

   function [blk,At,C,diagblkinfo,blockchange,parbarrier,X,Z] = ...
             detect_lblk(blk,At,C,b,parbarrier,X,Z,printlevel);

   if (nargin < 8); printlevel = 1; end
%%
%% Acum = abs(C) + sum_{k=1}^m abs(Ak)
%% but with diagonal elements removed. 
%%
   blkold = blk; 
   m = length(b);  
   ee = ones(m,1); 
   numdiagelt  = 0;  
   numblknew   = 0; 
   blockchange = zeros(size(blk,1),1);
   Acum        = cell(size(blk,1),1);    
   diagblkinfo = cell(size(blk,1),3); 
   for p=1:size(blk,1) 
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      if strcmp(pblk{1},'s') & (length(pblk{2}) == 1) & (size(At{p,1},2)==m)
         Acumtmp = smat(blk(p,:),abs(At{p})*ee,1) + abs(C{p});
         Acum{p} = Acumtmp - spdiags(diag(Acumtmp),0,n,n);          
         rownorm = sqrt(sum(Acum{p}.*Acum{p}))';
         idxdiag = find(rownorm < 1e-15);        
         if ~isempty(idxdiag)
            blockchange(p) = 1;
            numdiagelt = numdiagelt + length(idxdiag); 
            idxnondiag = setdiff([1:n]',idxdiag); 
            diagblkinfo{p,2} = idxdiag; 
            diagblkinfo{p,3} = idxnondiag; 
            if ~isempty(idxnondiag)
               numblknew = numblknew + 1; 
               diagblkinfo{p,1} = numblknew; 
            end
	 else
            numblknew = numblknew + 1; 
            diagblkinfo{p,1} = numblknew; 
         end
      else
         numblknew = numblknew + 1; 
         diagblkinfo{p,1} = numblknew; 
      end
   end
%%    
%% extract diagonal sub-blocks in nondiagonal blocks
%% into a single linear block
%%
   if any(blockchange)
      numblk = size(blkold,1); 
      idx_keepblk = []; 
      Atmp  = cell(1,m); 
      Adiag = cell(1,m); 
      C(numblk+1,1) = cell(1,1); 
      Cdiag = []; Xdiag = []; Zdiag = [];
      parbarrierdiag = [];
      for p = 1:numblk
         n = sum(blkold{p,2});
         if (blockchange(p)==1)
            idxdiag    = diagblkinfo{p,2}; 
            idxnondiag = diagblkinfo{p,3}; 
            if ~isempty(idxdiag); 
               blk{p,2} = length(idxnondiag);  
               len = length(idxdiag); 
               for k = 1:m;
                  Ak = mexsmat(blkold,At,1,p,k); 
                  tmp = diag(Ak); 
                  Atmp{k} = Ak(idxnondiag,idxnondiag);
                  Adiag{k} = [Adiag{k}; tmp(idxdiag)]; 
               end 
               tmp = diag(C{p,1}); 
               Cdiag = [Cdiag; tmp(idxdiag)];        
               C{p,1} = C{p,1}(idxnondiag,idxnondiag);
  	       At(p) = svec(blk(p,:),Atmp,1); 
               if (nargin >= 7)
                  parbarrierdiag = [parbarrierdiag, parbarrier{p}*ones(1,len)];
                  tmp = diag(X{p,1}); 
                  Xdiag = [Xdiag; tmp(idxdiag)]; 
                  tmp = diag(Z{p,1}); 
                  Zdiag = [Zdiag; tmp(idxdiag)];          
                  X{p,1} = X{p,1}(idxnondiag,idxnondiag);  
                  Z{p,1} = Z{p,1}(idxnondiag,idxnondiag);
               end
            end
            if ~isempty(idxnondiag)  
               idx_keepblk = [idx_keepblk, p];
	    else
               if (printlevel)
                  fprintf(' %2.0dth semidefinite block is actually diagonal\n',p); 
               end
            end
         else 
            idx_keepblk = [idx_keepblk, p];
         end
      end          
      blk{numblk+1,1} = 'l'; 
      blk{numblk+1,2} = numdiagelt;
      C{numblk+1,1}   = Cdiag;
      At(numblk+1,1) = svec(blk(numblk+1,:),Adiag,1); 
      idx_keepblk = [idx_keepblk, numblk+1];
      blk = blk(idx_keepblk,:); 
      C   = C(idx_keepblk,:); 
      At  = At(idx_keepblk,:);
      if (nargin >= 7)
         parbarrier{numblk+1,1} = parbarrierdiag;  
         X{numblk+1,1} = Xdiag; 
         Z{numblk+1,1} = Zdiag;
         parbarrier = parbarrier(idx_keepblk); 
         X = X(idx_keepblk); 
         Z = Z(idx_keepblk); 
      end
   end
%%******************************************************************



