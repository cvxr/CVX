%%*********************************************************
%% lmifun2: generate SDP data for the LMI constraint 
%%         of the form: 
%%
%%  [A*P*B'+B*P*A'  C*P*D']
%%  [D*P*C'          0    ]
%%
%%  P is the variable, P is symmetric.
%%
%%*********************************************************

   function Avec = lmifun2(A,B,C,D)

   n = size(A,2); n2 = n*(n+1)/2; 
   Acell = cell(1,n2); 
   m = size(D,1); 
   cnt = 1; 
%%
   ir2 = 1/sqrt(2); 
   for j = 1:n
      Bj = B(:,j); Aj = A(:,j); 
      Cj = C(:,j); Dj = D(:,j); 
      for i = 1:j 
         Ai = A(:,i); Bi = B(:,i); 
         Ci = C(:,i); Di = D(:,i); 
         if (i<j) 
            tmp = Ai*Bj' + Aj*Bi';
            tmp = tmp + tmp';
            tmp2 = Ci*Dj' + Cj*Di';  
            Acell{cnt} = ir2*sparse([tmp tmp2; tmp2' sparse(m,m)]); 
         else
            tmp = Ai*Bi'; 
            tmp = tmp + tmp'; 
            tmp2 = Ci*Di'; 
            Acell{cnt} = sparse([tmp tmp2; tmp2' sparse(m,m)]); 
         end
         cnt = cnt+1; 
      end
   end
%%
   blk{1,1} = 's'; blk{1,2} = size(A,1)+size(D,1); 
   Avec = svec(blk,Acell,1); 
%%*********************************************************
