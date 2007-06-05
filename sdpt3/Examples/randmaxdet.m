%%***************************************************************
%% randmaxdet: generate random determinant maximization 
%%             problem. 
%% min <C1,X1> + <C2,X2> - log(det(X2))
%%     A1(X1) + A2(X2) = b, X1 psd, X2 pd
%%
%% [blk,At,C,b,OPTIONS] = randmaxdet(n,p,m);
%% n = dimension of SDP variable X1
%% p = dimension of logdet variable X2
%% m = number of equality constraints
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%***************************************************************

  function [blk,At,C,b,OPTIONS] = randmaxdet(n,p,m);
   
   Y1 = randn(n);   Y1 = .5*(Y1 + Y1'); 
   Y1 = Y1 - min(0,1.1*min(real(eig(Y1))))*eye(n); 
   blk{1,1} = 's'; blk{1,2} = n;
   b = zeros(m,1); 
   F = cell(1,m); 
   for k = 1:m
      Fk = randn(n); Fk = .5*(Fk + Fk');
      F{1,k} = Fk; 
      b(k) = sum(sum(Y1.*Fk)); 
   end;
   At(1) = svec(blk(1,:),F);
   F0 = randn(n); F0 = .5*(F0+F0');
   C{1,1} = F0 - min(0,1.1*min(real(eig(F0))))*eye(n); 
   parbarrier{1} = 0; 
%%
   if (p > 0)
      Y2 = randn(p);   Y2 = .5*(Y2 + Y2'); 
      Y2 = Y2 - min(0,1.1*min(real(eig(Y2))))*eye(p);  
      blk{2,1} = 's'; blk{2,2} = p; 
      for k = 1:m
          Gk = randn(p); Gk = .5*(Gk + Gk');
          G{1,k} = Gk; 
          b(k) = b(k) + sum(sum(Y2.*Gk)); 
      end;
      At(2,1) = svec(blk(2,:),G);
      G0 = randn(p); G0 = .5*(G0+G0');
      C{2,1} = G0 - min(0,1.1*min(real(eig(G0))))*eye(p); 
      parbarrier{2} = 1; 
   end
   OPTIONS.parbarrier = parbarrier; 
%%***************************************************************
