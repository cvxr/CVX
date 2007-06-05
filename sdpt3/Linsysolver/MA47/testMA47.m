%%*********************************************************
%% testMA47: test whether MA47 has been installed 
%%           correctly. 
%%
%%
%%*********************************************************

  randn('seed',0);

  n = 100; 
  M = sprandn(n,n,0.2); M = M+M';
  rhs = ones(length(M),1); 

  [Lsymb,flag] = symbMA47(M); 
  [L,flag] = MA47fct(Lsymb,M);
  x = MA47slv(L,rhs); 

  xmatlab = M\rhs; 
  fprintf('\n norm(x-xmatlab)/norm(xmatlab) = %3.1e\n', ...
              norm(x-xmatlab)/norm(xmatlab)); 
%%*********************************************************
