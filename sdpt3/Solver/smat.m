%%*********************************************************
%% smat: compute the matrix smat(x).
%%
%%   M = smat(blk,x,isspM); 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**********************************************************

   function M = smat(blk,xvec,isspM);

   if (nargin < 3); isspM = zeros(size(blk,1),1); end 
%%
   if ~iscell(xvec) 
      if strcmp(blk{1},'s')
         M = mexsmat(blk,xvec,isspM);      
      else
         M = xvec; 
      end   
   else 
      M = cell(size(blk,1),1);     
      if (length(isspM)==1)
         isspM = isspM*ones(size(blk,1),1); 
      end
      for p=1:size(blk,1)
         pblk = blk(p,:);
         if strcmp(pblk{1},'s');
            M{p} = mexsmat(pblk,xvec{p},isspM(p));
         else
            M{p} = xvec{p}; 
         end   
      end
   end
%%*********************************************************

