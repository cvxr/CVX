%%***********************************************************
%% testchol: test file to verify that spchol
%%           has been set up correctly to solve
%%       
%%     X*y = b; 
%%
%%***********************************************************
 

  X = sprand(500,500,0.1); X = X*X'; X = 0.5*(X+X'); 
  m = length(X); 
  b = ones(m,1); b = sparse(b); %% create right-hand side vector
  fprintf('  testing sparse cholesky factorization rountines...\n'); 
%%
%% solve X*y = b
%%
  L = sparcholfun(symbcholfun(X),X);
  L.d(find(L.skip)) = inf;
  y = bwblkslvfun(L, fwblkslvfun(L,b) ./ L.d);
%%
%% verify results
%%
  res = norm(X(L.perm,L.perm)-L.L*spdiags(L.d,0,m,m)*L.L','fro'); 
  fprintf('  norm(X(perm,perm)- LDL^T) = %3.1e\n',res); 
  res = norm(X*y-b); 
  fprintf('  norm(X*y-b) = %3.1e\n',res); 
%%
%%***********************************************************
