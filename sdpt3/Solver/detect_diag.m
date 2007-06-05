%%*******************************************************************
%% detect_diag: detect diagonal blocks in the SDP data. 
%%
%% [blk,At,C,diagblkinfo] = detect_diag(blk,At,C,b); 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%******************************************************************

   function [blk,At,C,diagblkinfo] = detect_diag(blk,A,C,b,printlevel);

   if (nargin < 5); printlevel = 1; end

   m = length(b); 
   if iscell(A) & (size(A,2) == m)
      isvecA = 0;       
      At = svec(blk,A,1); 
   else
      isvecA = 1;
      At = A; 
   end
%%
%% concatenate diagonal blocks.
%%
   diagblk = [];  idx_keepblk = []; 
   for p = 1:size(blk,1)
      if strcmp(blk{p,1},'l'); 
         diagblk = [diagblk, p];  
      else
         idx_keepblk = [idx_keepblk, p]; 
      end 
   end
   if ~isempty(diagblk); 
      idx1 = diagblk(1); 
      idx_keepblk = [idx_keepblk, idx1];
   end
   if (length(diagblk) >= 2); 
      for p = 2:length(diagblk); 
         idxp = diagblk(p); 
   	 At{idx1} = [At{idx1}; At{idxp}]; 
         C{idx1,1} = [C{idx1,1}; C{idxp,1}]; 
         blk{idx1,2} = blk{idx1,2} + blk{idxp,2};             
      end
      blk = blk(idx_keepblk,:); 
      C = C(idx_keepblk,:); 
      At = At(idx_keepblk,:);
   end 
%%
%% Acum = abs(C) + sum_{k=1}^m abs(Ak)
%% but with diagonal elements removed. 
%%
   blkold = blk; 
   ee = ones(m,1); 
   numdiagelt  = 0;  
   blockchange = zeros(size(blk,1),1);
   Acum     = cell(size(blk,1),1);    
   diagblkinfo = cell(size(blk,1),1); 
   for p=1:size(blk,1) 
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      if strcmp(pblk{1},'s')
         Acumtmp = smat(blk(p,:),abs(At{p})*ee,1) + abs(C{p});
         Acum{p} = Acumtmp - spdiags(diag(Acumtmp),0,n,n);          
         rownorm = sum(Acum{p}.*Acum{p});
         diagblkinfo{p} = find(rownorm < 1e-13);
         if ~isempty(diagblkinfo{p})
            numdiagelt = numdiagelt + length(diagblkinfo{p}); 
            blockchange(p) = 1;
         end
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
      Cdiag = [];
      for p = 1:numblk
         n = sum(blkold{p,2});
         if strcmp(blkold{p,1},'s') 
            idxdiag    = diagblkinfo{p}; 
            idxnondiag = setdiff([1:n],idxdiag); 
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
               if (isvecA) 
  	          At(p) = svec(blk(p,:),Atmp,1);    
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
      C{numblk+1} = Cdiag;
      At(numblk+1,1) = svec(blk(numblk+1,:),Adiag,1); 
      idx_keepblk = [idx_keepblk, numblk+1];
      blk = blk(idx_keepblk,:); 
      C   = C(idx_keepblk,:); 
      At  = At(idx_keepblk,:);
   end
%%******************************************************************



