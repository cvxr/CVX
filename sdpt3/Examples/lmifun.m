%%*********************************************************
%% lmifun: generate SDP data for the LMI constraint 
%%         of the form: 
%%
%%  A*P*B' + B*P*A' + H*diag(d)*H' 
%%  P and d are variables, P is symmetric.
%%
%%  H is optional. 
%%*********************************************************

   function Avec = lmifun(A,B,H)

   n = size(A,2); n2 = n*(n+1)/2; r2 = sqrt(2); 
   Acell = cell(1,n2); 
   cnt = 1; 
%%
   for j = 1:n
      Bj = B(:,j); Aj = A(:,j); 
      for i = 1:j 
         Ai = A(:,i); Bi = B(:,i); 
         if (i<j) 
            tmp = Ai*Bj' + Aj*Bi'; 
            Acell{cnt} = sparse((tmp+tmp')/r2); 
         else
            tmp = Ai*Bi'; 
            Acell{cnt} = sparse(tmp+tmp'); 
         end
         cnt = cnt+1; 
      end
   end
%%
   if (nargin == 3) 
      for k = 1:size(H,2)
         Hk = H(:,k); 
         Acell{cnt} = sparse(Hk*Hk'); 
         cnt = cnt+1; 
      end    
   end
%%
   blktmp{1,1} = 's'; blktmp{1,2} = size(A,1); 
   Atmp = svec(blktmp,Acell,1); 
   Avec = Atmp{1};    
%%
%%*********************************************************
