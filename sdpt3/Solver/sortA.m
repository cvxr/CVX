%%*********************************************************************
%% sortA: sort columns of At{p} in ascending order according to the 
%%        number of nonzero elements. 
%%
%% [At,C,b,X0,Z0,permA,permZ] = sortA(blk,At,C,b,X0,Z0);
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*********************************************************************

   function [At,C,X0,Z0,permA,permZ] = sortA(blk,At,C,b,X0,Z0);

   global spdensity smallblkdim
%%
   if isempty(spdensity); spdensity = 0.4; end
   if isempty(smallblkdim); smallblkdim = 50; end
%%
   numblk = size(blk,1); 
   m = length(b); 
   nnzA = zeros(numblk,m); 
   permA = kron(ones(numblk,1),[1:m]); 
   permZ = cell(size(blk,1),1);
%%
   for p=1:size(blk,1)
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      numblk = length(pblk{2}); 
      if strcmp(pblk{1},'s') & (max(pblk{2}) > smallblkdim)
         n2 = sum(pblk{2}.*pblk{2});  n22 = sum(pblk{2}.*(pblk{2}+1))/2; 
         m1 = size(At{p,1},2);   
         if (length(pblk{2}) == 1)  
            tmp = abs(C{p}) + abs(Z0{p});
            if  (~isempty(At{p,1}))
                tmp = tmp + smat(blk(p,:),abs(At{p,1})*ones(m1,1),1);
            end
            if (nnz(tmp) < spdensity*n22); 
               per = symamd(tmp);    
               invper = zeros(n,1); invper(per) = [1:n]; 
               permZ{p} = invper;
               if (~isempty(At{p,1}))                  
                  isspAt = issparse(At{p,1});
                  for k = 1:m1
                     Ak = smat(pblk,At{p,1}(:,k),1); 
                     At{p,1}(:,k) = svec(pblk,Ak(per,per),isspAt); 
                  end
               end
               C{p}  = C{p}(per,per); 
               Z0{p} = Z0{p}(per,per); 
               X0{p} = X0{p}(per,per); 
            else
               per = [];
            end 
            if (length(pblk) > 2) & (~isempty(per)) 
               m2 = length(pblk{3}); 
               P = spconvert([(1:n)', per', ones(n,1)]);
               At{p,2} = P*At{p,2};
            end
         end
         if ~isempty(At{p,1}) & (mexnnz(At{p,1}) < m*n22/2)
            for k = 1:m1
                Ak = At{p,1}(:,k); 
                nnzA(p,k) = length(find(abs(Ak) > eps)); 
            end 
            [dummy,permAp] = sort(nnzA(p,1:m1)); 
            At{p,1}  = At{p,1}(:,permAp); 
            permA(p,1:m1) = permAp; 
         end
      elseif strcmp(pblk{1},'q') | strcmp(pblk{1},'l') | strcmp(pblk{1},'u'); 
         if ~issparse(At{p,1});
            At{p,1} = sparse(At{p,1}); 
         end
      end
   end
%%*********************************************************************
