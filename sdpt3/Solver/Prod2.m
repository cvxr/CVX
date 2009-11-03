%%*******************************************************************
%% Prod2: compute the block diagonal matrix A*B
%%
%%    C = Prod2(blk,A,B,options);
%%
%% INPUT:  blk = a cell array describing the block structure of A and B
%%         A,B = square matrices or column vectors.
%%
%%   options = 0  if no special structure 
%%             1  if C is symmetric
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

   function  C = Prod2(blk,A,B,options);
 
   global spdensity

   if (nargin == 3); options = 0; end;  
   iscellA = iscell(A); iscellB = iscell(B);
%%
   if (~iscellA & ~iscellB)
      if (size(blk,1) > 1); 
         error('Prod2: blk and A,B are not compatible'); 
      end;
      if strcmp(blk{1},'s')
         numblk = length(blk{2}); 
         isspA = issparse(A); isspB = issparse(B); 
         if (numblk > 1)
            if ~isspA; A=sparse(A); isspA=1; end
            if ~isspB; B=sparse(B); isspB=1; end
         end
         %%use_matlab = (options==0 & ~isspA & ~isspB) | (isspA & isspB); 
         use_matlab = (~isspA & ~isspB) | (isspA & isspB); 
         if (use_matlab)
            C = A*B;  
            if (options==1); C = 0.5*(C+C'); end; 
         else 
            C = mexProd2(blk,A,B,options);
         end
         checksparse = (numblk==1) & (isspA | isspB); 
         if (checksparse)
            n2 = sum(blk{2}.*blk{2}); 
            if (mexnnz(C) <= spdensity*n2);  
               if ~issparse(C); C = sparse(C); end; 
            else
               if issparse(C); C = full(C); end;
            end
         end              
      elseif (strcmp(blk{1},'q') | strcmp(blk{1},'l') | strcmp(blk{1},'u')) 
         C = A.*B;
      end 
   else
      error(['Prod2: A,B must be matrices']); 
   end  
%%*******************************************************************
