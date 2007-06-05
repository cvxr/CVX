%%*************************************************************
%% Asum: compute the matrix 
%%
%%  Ay = Asum(blk,A,y)
%%  
%%  input:  A = a CELL ARRAY with m columns. 
%%          y = mx1 vector.
%%          permAy = a permutation of [1:m] coding the order to 
%%                   sum the matrices y(k)*Ak, k = 1,...m. 
%%          iscmp = 1,  if Ay is complex
%%                  0,  otherwise.
%%                  
%%  output: Ay = sum_{k=1}^m y(k)*Ak, a column CELL ARRAY 
%%               with the same structure as A{:,1}. 
%%              
%% SDPT3: version 3.0
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 2 Feb 01
%%*************************************************************

   function Ay = Asum(blk,A,y)

   Ay = cell(size(blk,1),1); 
   m = length(y); 
   for p = 1:size(blk,1)
      n = sum(blk{p,2});  
      blktmp = blk{p,2};
      if strcmp(blk{p,1},'s'); 
         tmp = sparse(n,n); 
         for k = 1:m; tmp = tmp + A{p,k}*y(k); end; 
         if (length(blktmp) == 1)
            if (nnz(tmp) > 0.15*n*n);
               if issparse(tmp); tmp = full(tmp); end;
            else;
               if ~issparse(tmp); tmp = sparse(tmp); end;
            end;
         elseif (length(blktmp) > 1);
            if ~issparse(tmp); tmp = sparse(tmp); end;
         end;
      elseif strcmp(blk{p,1},'l'); 
         tmp = zeros(n,1); 
         for k = 1:m; tmp = tmp + A{p,k}*y(k); end; 
         if (nnz(tmp) > 0.15*n); 
            if issparse(tmp); tmp = full(tmp); end;
         else;
            if ~issparse(tmp); tmp = sparse(tmp); end;
         end;
      end;
      Ay{p} = tmp;
   end;
%%-------------------------------------------------------------



