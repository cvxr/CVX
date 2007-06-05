%%********************************************************************
%% scaling: scale the SDP data so that A_k,C,b have unit norm. 
%%
%%  [At,C,b,normA,normC,normb,X0,y0,Z0] = scaling(blk,At,C,b,X0,y0,Z0); 
%%
%%  The linear objective function value is unchanged under 
%%  the scaling. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%********************************************************************

   function [At,C,b,normA,normC,normb,X0,y0,Z0] = scaling(blk,At,C,b,X0,y0,Z0); 

   m = length(b); 
   numblk = size(blk,1);  
   
   normA = cell(numblk,1);  
   for p = 1:numblk
      pblk = blk(p,:);
      if strcmp(pblk{1},'s')
         m1 = size(At{p,1},2);
         if (m1 > 0)
            normAp2 = sum(sum(At{p,1}.*At{p,1})); 
	 else
	    normAp2 = 0;
         end
         if (length(pblk) > 2) %% for low rank constraints
            dd = At{p,3}; 
            m2 = m-m1; 
            ss = [0,cumsum(pblk{3})]; 
            for k=1:m2
               idx = [ss(k)+1:ss(k+1)]; 
               V = At{p,2}(:,idx); 
               ii = dd(idx,1)-ss(k); %% undo cumulative indexing
               jj = dd(idx,2)-ss(k); 
               len = pblk{3}(k); 
               D = spconvert([ii,jj,dd(idx,3); len,len,0]);
               tmp = V'*V*D; 
               normAp2 = normAp2 + sum(sum(tmp.*tmp')); 
            end
         end
      else
         normAp2 = sum(sum(At{p,1}.*At{p,1})); 
      end
      normAp = sqrt(normAp2);
      normA{p} = max(1,sqrt(normAp)); 
   end
%% 
   normb = max(1,norm(b));
   normC = 0; 
   for p = 1:numblk
      normC = max(normC,norm(C{p},'fro')); 
   end
   normC = max(1,normC); 
%%
   for p = 1:numblk
      pblk = blk(p,:);
      if strcmp(pblk{1},'s')
         m1 = size(At{p,1},2);
         m2 = m - m1;
         At{p,1} = At{p,1}/normA{p};
         if (m2 > 0) %% for low rank constaints
            At{p,3}(:,3) = At{p,3}(:,3)/normA{p};
         end
      else
         At{p,1} = At{p,1}/normA{p};
      end
      C{p}  = C{p}/(normC*normA{p}); 
      if (nargin == 7)
         X0{p} = X0{p}*normA{p}; 
         Z0{p} = Z0{p}/(normC*normA{p}); 
      end
   end   
   b = b/normb;
%%********************************************************************
