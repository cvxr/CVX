%%******************************************************************
%% dwd: distance weighted discrimination. 
%%
%% [w,beta,residp,residn,totalviolation,X,y,Z] = ...
%%             dwd(Ap,An,penalty);
%%
%% You're given two matrices, A+ and A-, whose columns
%% give points in R^n you want to separate "nicely" to allow 
%% classification of future points (e.g., to see if a patient 
%% has breast cancer or not). The "classical" solution puts a 
%% hyperplane between the points to maximize the minimum distance 
%% from a point to a hyperplane. Steve is interested in smaller 
%% instances, and is uncomfortable with the criterion, which is 
%% too subject to noise (like l_\infty fitting). He wants to use 
%% a criterion that depends on all the points. Initially he wanted 
%% to minimize the sum of the inverse squares of the distances, but 
%% he had no strong reason for the inverse square, so I suggested 
%% just the inverse distance. This turns into a nice socp! 
%% If there are mp points on the positive side and mn on the negative, 
%% the dimensions are: 2*(mp+mn)+1 constraints, 
%% (n+1) + (2) + (mp+mn)*3 variables in SOCs, with the dimensions as 
%% given (one large, the rest small), and mp+mn nonneg. variables. 
%% The latter are artificial variables to account for points that
%% are on the wrong side of the hyperplane: they can be moved a 
%% distance d at a cost penalty*d to put them on the right side.
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 20 Apr 02
%%******************************************************************

   function [w,beta,residp,residn,totalviolation,X,y,Z] = ...
             dwd(Ap,An,penalty);

   [np,mp] = size(Ap);
   [nn,mn] = size(An);
   n = np;
   nv = 1 + n + 2 + 3*(mp + mn) + mp + mn;
   nc = 1 + 2*(mp + mn);
%%
   blk = cell(2,2);
   blk{1,1} = 'q';
   blk{1,2} = [n+1, 2, 3*ones(1, mp + mn)];
   blk{2,1} = 'l';
   blk{2,2} = mp + mn;
%%
   Avec = cell(2,1);
   A = sparse(nc,nv-mp-mn);
   A(1:mp,2:n+3) = [Ap',zeros(mp,1),ones(mp,1)];
   A(mp+1:mp+mn,2:n+3) = -[An',zeros(mn,1),ones(mn,1)];
   A(1:mp,n+4:3:n+3+3*mp) = - speye(mp,mp);
   A(1:mp,n+6:3:n+5+3*mp) = - speye(mp,mp);
   A(mp+1:mp+mn,3*mp+n+4:3:3*mp+n+3+3*mn) = - speye(mn,mn);
   A(mp+1:mp+mn,3*mp+n+6:3:3*mp+n+5+3*mn) = - speye(mn,mn);
   A(mp+mn+1,1) = 1;
   A(mp+mn+2:mp+mn+1+mp,n+5:3:n+4+3*mp) = speye(mp,mp);
   A(mp+mn+1+mp+1:mp+mn+1+mp+mn,3*mp+n+5:3:3*mp+n+4+3*mn) = speye(mn,mn);
%%   
   Avec{1,1} = A;
   Avec{2,1} = [speye(mp+mn,mp+mn);zeros(1+mp+mn,mp+mn)];
   b = [zeros(mp+mn,1);ones(1+mp+mn,1)];
%%
   C = cell(2,1);
   c = zeros(nv-mp-mn,1);
   c(n+4:3:n+3+3*mp) = ones(mp,1);
   c(n+6:3:n+5+3*mp) = -ones(mp,1);
   c(3*mp+n+4:3:3*mp+n+3+3*mn) = ones(mn,1);
   c(3*mp+n+6:3:3*mp+n+5+3*mn) = -ones(mn,1);
   C{1,1} = c;
   C{2,1} = penalty*ones(mp+mn,1);
%%
   [obj,X,y,Z] = sqlp(blk,Avec,C,b);
   X1 = X{1};
   omega = X1(1);
   w = X1(2:n+1);
   alpha = X1(n+2);
   beta  = X1(n+3);
   residp = Ap'*w + beta*ones(mp,1);
   residn = An'*w + beta*ones(mn,1);
   X2 = X{2};
   totalviolation = sum(X2);
%%******************************************************************

