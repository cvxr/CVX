%%******************************************************
%% MA47fct: factorize a symmetric indefinite matrix
%%          given the symbolic factor Lsymb. 
%%
%%******************************************************

   function  [L,flag] = MA47fct(Lsymb,X);

   [L.a,L.iw,flag] = mexMA47fct(X,Lsymb,[1e-10,0]); 
   L.perm = [1:length(X)]; 
%%******************************************************
