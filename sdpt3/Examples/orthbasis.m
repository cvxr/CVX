%%****************************************************************
%% orthbasis: Compute an orthonomal basis from 
%%            the power basis {I,A,A^2, .... ,A^m}
%%            via an Arnoldi-type iteration. 
%%
%%   [Q,H,C] = orthbasis(A,m);
%% 
%%   Output:  Q = a cell array containing the orthonormal basis
%%                obtained from {I,A,A^2, .... ,A^m}. 
%%
%% use in chebymat.m, igmres.m
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%****************************************************************

  function [Q,H,C] = orthbasis(A,m);

  N = length(A); 
  C(1,1) = 1/sqrt(N); 
  Q = cell(1,m+1); 
  Q{1} = eye(N)/sqrt(N); 
  for k = 1:m
      V = A*Q{k}; 
      Vold = V; 
      for j = 1:k
          H(j,k) = sum(sum(conj(Q{j}).*V)); 
          V = V - H(j,k)*Q{j};
      end; 
      if (norm(V,'fro') < norm(Vold,'fro')); 
         for j = 1:k
            s(j,1) = sum(sum(conj(Q{j}).*V));
            V = V - s(j)*Q{j};   
         end; 
         H(1:k,k) = H(1:k,k) + s; 
      end; 
      nrm = norm(V,'fro'); 
      H(k+1,k) = nrm;
      Q{k+1} = V/nrm; 
      C(1:k+1,k+1) = ([0; C(1:k,k)] - [C(1:k,1:k)*H(1:k,k); 0])/nrm;  
  end; 
%%================================================================



