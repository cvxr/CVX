%%**************************************************
%% symbMA47: compute symbolic factorization of
%%           a symmetric indefinite matrix.
%%
%%**************************************************

  function [Lsymb,flag] = symbMA47(X)

  nnztmp = (mexnnz(X) + mexnnz(diag(X)))/2; 
  [keep,Jcn,flag,info] = mexMA47syb(nnztmp,X);
  if (flag)
     Lsymb = []; 
  else
     Lsymb.keep = keep;    Lsymb.Jcn = Jcn;   
     Lsymb.La = 2*info(1); Lsymb.Liw = 2*info(2);  
  end
%%**************************************************
