%%*****************************************************
%% viewgraph: plot the graph of an adjacency matrix.
%%
%% SDPT3: version 3.0 
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last modified: 2 Feb 01
%%*****************************************************

   function viewgraph(B);

   n = length(B); 
   z = exp(2*pi*sqrt(-1)*[0:n-1]/n);
   circle = exp(2*pi*sqrt(-1)*linspace(0,1,100)); 
   plot(circle,':'); hold on
   axis('square'); 
   for i = 1:n
       idx = find(B(i,i+1:n) == 1);
       idx = idx + i;
       for k = 1:length(idx)
           z1 = z(i); z2 = z(idx(k)); 
           plot([real(z1) real(z2)],[imag(z1) imag(z2)]); 
       end;
   end;
   plot(z,'.r','markersize',16); 
   hold off;
%%=====================================================
