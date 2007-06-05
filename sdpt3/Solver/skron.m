%%***********************************************************************
%% skron:  Find the matrix presentation of
%%         symmetric kronecker product skron(A,B), where
%%         A,B are symmetric.
%%
%% Important: A,B are assumed to be symmetric. 
%%
%% K = skron(blk,A,B); 
%%
%% blk: a cell array specifying the block diagonal structure of A,B. 
%% 
%% (ij)-column of K = 0.5*svec(AUB + BUA),  where
%%                     U  = xij*(ei*ej' + ej*ei')
%%                    xij = 1/2       if i=j
%%                        = 1/sqrt(2) otherwise.
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%***********************************************************************
 
  function K = skron(blk,A,B); 

  if iscell(A) & ~ iscell(B)
     error('skron: A,B must be both matrices or both cell arrays')
  end
  if iscell(A)  
     K = cell(size(blk,1),1); 
     for p = 1:size(blk,1) 
        if (norm(A{p}-B{p},'fro') < 1e-13)
           sym = 1; 
        else
           sym = 0; 
        end     
        if strcmp(blk{p,1},'s')
           K{p} = mexskron(blk(p,:),A{p},B{p},sym);
        end
     end
  else
     if (norm(A-B,'fro') < 1e-13)
        sym = 1; 
     else
        sym = 0; 
     end     
     if strcmp(blk{1,1},'s')
        K = mexskron(blk,A,B,sym);
     end
  end
%%***********************************************************************
