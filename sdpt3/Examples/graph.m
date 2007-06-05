%%******************************************************
%% graph: generate random adjacency matrix.
%%      
%%  B = graph(n,prob);
%%
%%  see maxcut.m 
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%******************************************************

   function B = graph(n,prob);

   B = zeros(n,n);
   if nargin <= 1; prob = 0.5; end;
   for i = 1:n
       for j = i+1:n 
           r = rand(1); 
           if (r < prob); B(i,j) = 1; B(j,i) = 1;  end;
       end;
   end;
%%======================================================
