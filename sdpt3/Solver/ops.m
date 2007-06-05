%%******************************************************************
%% ops:  
%%
%%   Z = ops(X,operand,Y,alpha); 
%%
%%  INPUT:        X = a matrix or a scalar
%%                    or a CELL ARRAY consisting only of matrices 
%%          operand = sym, transpose, triu, tril,
%%                    real, imag, sqrt, abs, max, min, nnz,
%%                    spdiags, ones, norm, sum, row-norm, 
%%                    rank1, rank1inv, inv
%%                    +,  -, *, .*,  ./, .^ 
%%     Y (optional) = a matrix or a scalar 
%%                    or a CELL ARRAY consisting only of matrices
%% alpha (optional) = a scalar
%%                    or the variable blk. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%******************************************************************

   function Z = ops(X,operand,Y,alpha); 

   global spdensity 

   if nargin == 2; 
      if strcmp(operand,'sym'); 
         if ~iscell(X); 
            [m,n] = size(X); 
            if (m == n); 
               Z = (X+X')/2; 
            elseif (n == 1);  
               Z = X;
            else; 
               error('X must be square matrix or a column vector'); 
            end; 
         else
            Z = cell(size(X)); 
            for p = 1:length(X); 
                [m,n] = size(X{p}); 
                if (m == n); 
                   Z{p} = (X{p}+X{p}')/2; 
                elseif (n == 1); 
                   Z{p} = X{p};
                else; 
                   error('X{p} must be square matrix or a column vector'); 
                end; 
            end;
         end;
      elseif strcmp(operand,'sqrt') | strcmp(operand,'abs') | ...
         strcmp(operand,'real') | strcmp(operand,'imag');
         if ~iscell(X); 
            eval(['Z = ',operand,'(X);']);  
         else;
            Z = cell(size(X)); 
            for p = 1:length(X); 
               eval(['Z{p} = ',operand,'(X{p});']);  
            end;
         end;
      elseif strcmp(operand,'max') | strcmp(operand,'min') | ...
             strcmp(operand,'sum'); 
         if ~iscell(X); 
            eval(['Z = ',operand,'(X);']);  
         else;
            Z = []; 
            for p = 1:length(X); 
                eval(['Z = [Z  ',operand,'(X{p})','  ];']);  
            end;
         end; 
         eval(['Z = ',operand,'(Z);']);  
      elseif strcmp(operand,'transpose') | strcmp(operand,'triu') | ...
             strcmp(operand,'tril'); 
         if ~iscell(X); 
            if (size(X,1) == size(X,2)); 
               eval(['Z = ',operand,'(X);']);
            elseif (size(X,2) == 1); 
               eval(['Z = X;']);
            else; 
               error('X must be square matrix or a column vector'); 
            end; 
         else
            Z = cell(size(X)); 
            for p = 1:length(X); 
                if (size(X{p},1) == size(X{p},2)); 
                   eval(['Z{p} = ',operand,'(X{p});']); 
                elseif (size(X{p},2) == 1); 
                   eval(['Z{p} = X{p};']);
                else; 
                   error('X{p} must be square matrix or a column vector'); 
                end; 
            end;
         end;
      elseif strcmp(operand,'norm');
         if ~iscell(X); 
            Z = full(sqrt(sum(sum(X.*X))));
         else
            Z = 0; 
            for p = 1:length(X); Z = Z + sum(sum(X{p}.*X{p})); end;
            Z = sqrt(Z); 
         end;
      elseif strcmp(operand,'inv'); 
         if ~iscell(X);
            [m,n] = size(X); n2 = n*n; 
            Z = inv(X); 
            if (nnz(Z) > spdensity*n2) & issparse(Z); 
               Z = full(Z); 
            elseif (m > 1 & n == 1);
               Z = 1./X; 
               if (nnz(Z) > spdenssity*n) & issparse(Z); 
                  Z = full(Z); 
               end;
            end;
         else
            Z = cell(size(X)); 
            for p = 1:length(X); 
               [m,n] = size(X{p}); n2 = n*n;
               if (m == n); 
                  Z{p} = inv(X{p}); 
                  if (nnz(Z{p}) > spdensity*n2) & issparse(Z{p}); 
                     Z{p} = full(Z{p});  
                  end
               elseif (m > 1 & n == 1);
                  Z{p} = 1./X{p}; 
                  if (nnz(Z{p}) > spdensity*n) & issparse(Z{p}); 
                     Z{p} = full(Z{p}); 
                  end;
               end;
            end; 
         end;
      elseif strcmp(operand,'getM'); 
         if ~iscell(X); 
            Z = size(X,1);
         else
            for p = 1:length(X); Z(p) = size(X{p},1); end;
            Z = sum(Z); 
         end; 
      elseif strcmp(operand,'nnz');
         if ~iscell(X); 
            Z = nnz(X);  
         else;
            for p = 1:length(X); 
               Z(p) = nnz(X{p}); 
            end;
            Z = sum(Z);  
         end; 
      elseif strcmp(operand,'ones');
         if ~iscell(X); 
            Z = ones(size(X));
         else; 
            Z = cell(size(X)); 
            for p = 1:length(X);
                Z{p} = ones(size(X{p}));
            end;
         end;
      elseif strcmp(operand,'row-norm'); 
         if ~iscell(X);
            if (size(X,2) == size(X,1)); 
               Z = sqrt(sum((X.*conj(X))'))';
            elseif (size(X,2) == 1);  
               Z = abs(X); 
            end;
         else; 
            Z = cell(size(X)); 
            for p = 1:length(X);
               if (size(X{p},2) == size(X{p},1)); 
                  Z{p} = sqrt(sum((X{p}.*conj(X{p}))'))'; 
               elseif (size(X{p},2) == 1);  
                  Z{p} = abs(X{p}); 
               end;
            end;
         end;        
      end; 
   end;
%%
   if (nargin == 3)
      if strcmp(operand,'spdiags');
         if ~iscell(Y); 
            [m,n] = size(Y); 
            if (m == n); 
               Z = spdiags(X,0,m,n);
            else
               Z = X;
            end
         else 
            Z = cell(size(Y)); 
            for p = 1:length(Y);
               [m,n] = size(Y{p}); 
               if (m == n); 
                  Z{p} = spdiags(X{p},0,m,n);
               else;
                  Z{p} = X{p};
               end
            end
         end
      elseif strcmp(operand,'inprod')
         if ~iscell(X) & ~iscell(Y)
    	    Z = (Y'*X)';
	 elseif iscell(X) & iscell(Y)
    	    Z = zeros(size(X{1},2),1); 
   	    for p=1:length(X)
	        Z = Z + (Y{p}'*X{p})'; 
            end
         end
      elseif strcmp(operand,'+') | strcmp(operand,'-') | ...
             strcmp(operand,'/') | strcmp(operand,'./') | ...
             strcmp(operand,'*') | strcmp(operand,'.*') | ...
             strcmp(operand,'.^');
         if (~iscell(X) & ~iscell(Y)); 
            eval(['Z = X',operand,'Y;']); 
         elseif (iscell(X) & iscell(Y))
            Z = cell(size(X)); 
            for p = 1:length(X); 
 	       if (size(X{p},2) == 1) & (size(Y{p},2) == 1) & ... 
                  (strcmp(operand,'*') | strcmp(operand,'/')); 
                  eval(['Z{p} = X{p}.',operand,'Y{p};']); 
               else
                  eval(['Z{p} = X{p} ',operand,'Y{p};']);                   
               end
            end 
         elseif (iscell(X) & ~iscell(Y)); 
            Z = cell(size(X)); 
            for p = 1:length(X); 
                eval(['Z{p} = X{p}',operand,'Y;']); 
            end
         elseif (~iscell(X) & iscell(Y)); 
            Z = cell(size(Y)); 
            for p = 1:length(Y); 
                eval(['Z{p} = X',operand,'Y{p};']); 
            end
         end
      else
         error([operand,' is not available, check input arguments']); 
      end
   end
%%
   if (nargin == 4)
      if strcmp(operand,'rank1') | strcmp(operand,'rank1inv'); 
         Z = cell(size(alpha,1),1); 
         for p = 1:size(alpha,1);
            if ~strcmp(alpha{p,1},'diag');
               blktmp = alpha{p,2}; 
               if (length(blktmp) == 1); 
                  if strcmp(operand,'rank1');                
                     Z{p} = (X{p}*Y{p}' + Y{p}*X{p}')/2; 
                  else;
                     Z{p} = 2./(X{p}*Y{p}' + Y{p}*X{p}'); 
                  end
               else             
                  Xp = X{p}; 
                  Yp = Y{p};   
                  n = sum(blktmp);             
                  Zp = sparse(n,n);
                  s = [0 cumsum(blktmp)];  
                  if strcmp(operand,'rank1');
                     for i = 1:length(blktmp)
                        pos = [s(i)+1 : s(i+1)]; 
                        x = Xp(pos); 
                        y = Yp(pos); 
                        Zp(pos,pos) = sparse((x*y' + y*x')/2);
                     end;
                     Z{p} = Zp; 
                  else 
                     for i = 1:length(blktmp)
                        pos = [s(i)+1 : s(i+1)]; 
                        x = Xp(pos); 
                        y = Yp(pos); 
                        Zp(pos,pos) = sparse(2./(x*y' + y*x'));
                     end
                     Z{p} = Zp; 
                  end
               end
            elseif strcmp(alpha{p,1},'diag'); 
               if strcmp(operand,'rank1'); 
                  Z{p} = X{p}.*Y{p};
               else
                  Z{p} = 1./(X{p}.*Y{p});
               end
            end
         end
      elseif strcmp(operand,'+') | strcmp(operand,'-');
         if ~iscell(X) & ~iscell(Y); 
            eval(['Z = X',operand,'alpha*Y;']); 
         elseif (iscell(X) & iscell(Y)); 
            Z = cell(size(X)); 
            for p = 1:length(X); 
               eval(['Z{p} = X{p}',operand,'alpha*Y{p};']); 
            end
         else
            error('X, Y are different objects'); 
         end
      else
         error([operand,' is not available']); 
      end
   end
%%============================================================



