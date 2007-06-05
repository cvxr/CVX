%*********************************************************
%% svec: compute the vector svec(M),
%%
%%   x = svec(blk,M,isspx);
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**********************************************************

  function x = svec(blk,M,isspx); 
     
   if iscell(M) 
      if (size(blk,1) ~= size(M,1))
         error('svec: number of rows in blk and M not equal');
      end
      if (nargin == 2)
         %%if (size(M,2) == 1)
         %%   isspx = zeros(size(blk,1),1); 
         %%else 
         %%   isspx = ones(size(blk,1),1); 
         %%end
         isspx = ones(size(blk,1),1); 
      else
         if (length(isspx) < size(blk,1))
            isspx = ones(size(blk,1),1); 
         end
      end
      x = cell(size(blk,1),1); 
      for p=1:size(blk,1)
         pblk = blk(p,:);  
         n = sum(pblk{2});  m = size(M,2); 
         if strcmp(pblk{1},'s')
            n2 = sum(pblk{2}.*(pblk{2}+1))/2; 
            if (isspx(p)); 
               x{p} = sparse(n2,m); 
            else 
               x{p} = zeros(n2,m); 
            end
            numblk = length(pblk{2}); 
            if (pblk{2} > 0)
               for k = 1:m
                  if (numblk > 1) & ~issparse(M{p,k});
                     x{p}(:,k) = mexsvec(pblk,sparse(M{p,k}),isspx(p)); 
                  else
                     x{p}(:,k) = mexsvec(pblk,M{p,k},isspx(p)); 
                  end
               end               
            end
	 else
            if (isspx(p)) 
               x{p} = sparse(n,m); 
            else 
               x{p} = zeros(n,m); 
            end
            for k = 1:m 
                x{p}(:,k) = M{p,k};
            end
         end
      end
   else 
      if strcmp(blk{1},'s')
         numblk = length(blk{2}); 
         if (numblk > 1) & ~issparse(M);
            x = mexsvec(blk,sparse(M),1); 
         else
            x = mexsvec(blk,sparse(M)); 
         end
      else
         x = M;
      end
   end 
%%**********************************************************   

