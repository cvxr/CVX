%%***************************************************************************
%% blkeig: compute eigenvalue decomposition of a cell array
%%         whose contents are square matrices or the diagonal
%%         of a diagonal matrix. 
%% 
%% [d,V] = blkeig(blk,X);
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%***************************************************************************

    function [d,V] = blkeig(blk,X);

    spdensity = 0.5; 

    if ~iscell(X); 
       if strcmp(blk{1},'s'); 
          blktmp = blk{2}; 
          if (length(blktmp) == 1);
             if (nargout == 1); 
                d = eig(full(X)); 
             elseif (nargout == 2);
                [V,d] = eig(full(X));
                d = diag(d);  
                [d,idx] = sort(d);
                V = V(:,idx);  
             end
          else
             if (nargout == 2); 
                V = sparse(length(X),length(X));  
             end;
             s = [0 cumsum(blktmp)]; 
             for i = 1:length(blktmp)
                pos = [s(i)+1 : s(i+1)];   
                if (nargout == 1); 
                   lam = eig(full(X(pos,pos))); 
                elseif (nargout == 2);           
                   [evec,lam] = eig(full(X(pos,pos))); 
                   lam  = diag(lam);
                   V(pos,pos) = sparse(evec); 
                end; 
                d(pos,1) = lam;
             end
          end
          n2 = sum(blktmp.*blktmp); 
          if (nargout == 2); 
             if (nnz(V) <= spdensity*n2);
                V = sparse(V);  
             else
                V = full(V); 
             end; 
          end; 
       elseif strcmp(blk{1},'l');
          if (nargout == 2); 
             V = ones(size(X)); d = X; 
          elseif (nargout == 1); 
             d = X; 
          end
       end
    else
       if (nargout == 2); 
          V = cell(size(X));  d = cell(size(X));  
          for p = 1:size(blk,1); 
	     [d{p},V{p}] = blkeig(blk(p,:),X{p}); 
          end;
       elseif (nargout == 1); 
          d = cell(size(X));  
          for p = 1:size(blk,1); 
	     d{p} = blkeig(blk(p,:),X{p}); 
          end
       end
    end
%%***************************************************************************
