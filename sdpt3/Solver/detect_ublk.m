%%*******************************************************************
%% detect_ublk: search for implied free variables in linear
%%              block. 
%% [blk2,At2,C2,ublkinfo] = detect_ublk(blk,At,C); 
%%
%% i1,i2: indices corresponding to splitting of unrestricted varaibles
%% i3   : remaining indices in the linear block
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

   function [blk2,At2,C2,ublkinfo,parbarrier2,X2,Z2] = ...
             detect_ublk(blk,At,C,parbarrier,X,Z,printlevel); 
   
   if (nargin < 7); printlevel = 1; end
  
   blk2 = blk; At2 = At; C2 = C; 
   if (nargin >= 6)
      parbarrier2 = parbarrier; 
      X2 = X; Z2 = Z; 
   else
      X2 = []; Z2 = [];
   end
   numblk = size(blk,1);
   ublkinfo = cell(size(blk,1),3);  
   tol = 1e-14;
%%
   numblknew = numblk; 
%%
   for p = 1:numblk
      pblk = blk(p,:);
      m = size(At{p},2);        
      if strcmp(pblk{1},'l')
         r = randmat(1,m,0,'n');
         stime = cputime;
         Ap = At{p}'; Cp = C{p};
	 ApTr = (r*Ap)';
	 [sApTr,perm] = sort(abs(ApTr));
	 idx0 = find(abs(diff(sApTr)) < tol);
         i1 = []; i2 = [];
	 if ~isempty(idx0)
            n = pblk{2}; 
            i1 = perm(idx0); i2 = perm(idx0+1);
	    Api1 = Ap(:,i1);
	    Api2 = Ap(:,i2);
	    Cpi1 = Cp(i1)';
	    Cpi2 = Cp(i2)';
	    idxzr = find(abs(Cpi1+Cpi2) < tol & sum(abs(Api1+Api2),1) < tol);
            if ~isempty(idxzr)
               i1 = i1(idxzr');
               i2 = i2(idxzr');
               blk2{p,1} = 'u'; 
               blk2{p,2} = length(i1); 
               At2{p} = Ap(:,i1)'; 
               C2{p}  = Cp(i1); 
               if (printlevel)
                   fprintf('\n %1.0d linear variables from unrestricted variable.\n',...
                            2*length(i1)); 
               end
               if (nargin >= 6)
                  parbarrier2{p} = parbarrier{p}(i1); 
                  X2{p} = X{p}(i1)-X{p}(i2); 
                  Z2{p} = zeros(length(i1),1); 
               end
               i3 = setdiff([1:n],union(i1,i2));
               if ~isempty(i3)
                  numblknew = numblknew + 1; 
                  blk2{numblknew,1} = 'l'; 
                  blk2{numblknew,2} = length(i3); 
                  At2{numblknew,1}  = Ap(:,i3)'; 
                  C2{numblknew,1}   = Cp(i3); 
                  if (nargin >= 6)
                     parbarrier2{numblknew,1} = parbarrier{p}(i3); 
                     X2{numblknew,1} = X{p}(i3); Z2{numblknew,1} = Z{p}(i3); 
                  end               
               end
	       ublkinfo{p,1} = i1; ublkinfo{p,2} = i2; ublkinfo{p,3} = i3;
            end
         end
      end
   end
%%*******************************************************************
