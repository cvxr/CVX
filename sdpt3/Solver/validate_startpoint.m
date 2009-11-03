%%***********************************************************************
%% validate_startpoint: validate_startpoint starting point X0,y0,Z0
%%
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%***********************************************************************

   function [X0,Z0] = validate_startpoint_startpoint(blk,X0,Z0,spdensity,iscmp); 

   if (nargin < 5); iscmp = 0; end
   if (nargin < 4); spdensity = 0.4; end   
%%
   if ~iscell(X0) | ~iscell(Z0); 
      error('validate_startpoint: X0, Z0 must be cell arrays'); 
   end
   if (min(size(X0))~=1 | min(size(Z0))~=1); 
      error('validate_startpoint: cell array X, Z can only have 1 column or row'); 
   end
   if (size(X0,2) > size(X0,1)); X0 = X0'; end;
   if (size(Z0,2) > size(Z0,1)); Z0 = Z0'; end;
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      n2 = sum(pblk{2}.*pblk{2});  
      numblk = length(pblk{2});
      if strcmp(pblk{1},'s');
         if (iscmp)
            X0{p} = [real(X0{p}),-imag(X0{p}); imag(X0{p}), real(X0{p})]; 
            Z0{p} = [real(Z0{p}),-imag(Z0{p}); imag(Z0{p}), real(Z0{p})]; 
         end
         if ~all(size(X0{p}) == n) | ~all(size(Z0{p}) == n); 
            error('validate_startpoint: blk and X0,Z0 are not compatible'); 
         end 
         if (norm([X0{p}-X0{p}' Z0{p}-Z0{p}'],inf) > 2e-13); 
            error('validate_startpoint: X0,Z0 not symmetric'); 
         end
         if (nnz(X0{p}) < spdensity*n2) | (numblk > 1) ; 
            if ~issparse(X0{p}); X0{p} = sparse(X0{p}); end; 
         else
            if issparse(X0{p}); X0{p} = full(X0{p}); end;
         end
         if (nnz(Z0{p}) < spdensity*n2) | (numblk > 1); 
            if ~issparse(Z0{p}); Z0{p} = sparse(Z0{p}); end; 
         else
            if issparse(Z0{p}); Z0{p} = full(Z0{p}); end;
         end
      elseif strcmp(pblk{1},'q') | strcmp(pblk{1},'l') | strcmp(pblk{1},'u'); 
         if ~all([size(X0{p},2) size(Z0{p},2)]==1); 
            error(['validate_startpoint: ',num2str(p),...
            '-th block of X0,Z0 must be column vectors']);
         end
         if ~all([size(X0{p},1) size(Z0{p},1)]==n); 
            error(['validate_startpoint: blk, and X0,Z0, are not compatible']); 
         end              
         if (nnz(X0{p}) < spdensity*n); 
            if ~issparse(X0{p}); X0{p} = sparse(X0{p}); end; 
         else
            if issparse(X0{p}); X0{p} = full(X0{p}); end;
         end
         if (nnz(Z0{p}) < spdensity*n); 
            if ~issparse(Z0{p}); Z0{p} = sparse(Z0{p}); end; 
         else
            if issparse(Z0{p}); Z0{p} = full(Z0{p}); end;
         end
	 if strcmp(pblk{1},'l') & (any(X0{p} < 1e-12) | any(Z0{p} < 1e-12))
            error(['X0 or Z0 is not in nonnegative cone']); 
         end
         if strcmp(pblk{1},'q'); 
            s = 1+[0, cumsum(pblk{2})]; len = length(pblk{2}); 
            if any(X0{p}(s(1:len)) < 1e-12) | any(Z0{p}(s(1:len)) < 1e-12)
               error(['X0 or Z0 is not in socp cone']); 
            end
         end
      end
   end 
%%***********************************************************************
