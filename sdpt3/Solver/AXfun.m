%%*************************************************************************
%% AXfun: compute AX(k) = <Ak,X>, k = 1:m
%%
%%   AX = AXfun(blk,At,permA,X);
%%
%% Note: permA may be set to [] if no permutation is neccessary. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*************************************************************************

  function AX = AXfun(blk,At,permA,X);
  
  if isempty(permA); ismtpermA = 1; else; ismtpermA = 0; end

  for p = 1:size(blk,1); 
     pblk = blk(p,:);
     if strcmp(pblk{1},'s')
        m1 = size(At{p,1},2); 
        if (p==1) 
           if (length(pblk) > 2); m2 = length(pblk{3}); else; m2 = 0; end
           m = m1 + m2;
           AX = zeros(m,1);  tmp = zeros(m,1);
        end
        if (~isempty(At{p,1}))
           if (ismtpermA)
              tmp = (svec(pblk,X{p})'*At{p,1})';
              %%tmp = mexinprod(blk,At,svec(pblk,X{p}),m1,p); 
           else
              tmp(permA(p,1:m1),1) = (svec(pblk,X{p})'*At{p,1})'; 
              %%tmp(permA(p,1:m1),1) = mexinprod(blk,At,svec(pblk,X{p}),m1,p);
           end
        end
        if (length(pblk) > 2)  %% for low rank constraints
	   m2 = length(pblk{3});
           dd = At{p,3};
           len = sum(pblk{3});
           DD = spconvert([dd(:,2:4); len,len,0]);    
           XVD = X{p}*At{p,2}*DD;
	   if (length(X{p}) > 1)
              tmp2 = sum(At{p,2}.*XVD)'; 
           else
              tmp2 = (At{p,2}.*XVD)'; 
           end
           tmp(m1+[1:m2]) = mexqops(pblk{3},tmp2,ones(length(tmp2),1),1); 
        end
        AX = AX + tmp;
     else
        if (p==1); m = size(At{p,1},2); AX = zeros(m,1);  tmp = zeros(m,1); end
        AX = AX + (X{p}'*At{p,1})'; 
     end
  end
%%*************************************************************************
