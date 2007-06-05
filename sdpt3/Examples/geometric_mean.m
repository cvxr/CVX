%%******************************************************************
%% geometric_mean: an example with rotated cone variables
%%
%%  max {prod(d+B*x) : d+B*x > 0, x <= 10}  
%% 
%%  where B = 4xn matrix, 
%%        d = 4x1 vector
%%
%% [blk,At,C,b,xx] = geometric_mean(B,d,solve); 
%%
%% E.g. p = 6; m = 4; B = [rand(2,p); -rand(2,p)]; d = rand(m,1);
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%******************************************************************

   function [blk,At,C,b,xx] = geometric_mean(B,d,solve); 

   if (nargin == 2); solve = 0; end

   n = size(B,2); 
   zz = zeros(1,n); 
   r2 = sqrt(2); 
   blk{1,1} = 'r'; blk{1,2} = [3,3,3];
   At{1,1} = -[B(1,:),0,0,0; B(2,:),0,0,0; zz,r2,0,0; ...
               B(3,:),0,0,0; B(4,:),0,0,0; zz,0,r2,0; ...
	       zz, 1,0,0; zz, 0,1,0; zz, 0,0,r2]; 
   C{1,1} = [d(1:2); 0; d(3:4); 0; 0;0;0]; 
   b = [zeros(n,1); 0;0;1]; 
   blk{2,1} = 'l'; blk{2,2} = n; 
   At{2,1} = [eye(n), zeros(n,3)]; C{2,1} = 10*ones(n,1);
   if (solve)
      [bblk,AAt,CC,bb,T] = convertRcone(blk,At,C,b);
      [obj,X,y,Z] = sqlp(bblk,AAt,CC,bb);
      xx = y(1:n); 
   end
%%******************************************************************
