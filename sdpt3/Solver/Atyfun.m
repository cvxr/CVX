%%*********************************************************
%% Atyfun: compute sum_{k=1}^m yk*Ak. 
%%
%%  Q = Atyfun(blk,At,permA,isspAy,y);
%%
%% Note: permA and isspAy may be set to [].
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**********************************************************

  function Q = Atyfun(blk,At,permA,isspAy,y);

  if isempty(permA); ismtpermA = 1; else; ismtpermA = 0; end

  Q = cell(size(blk,1),1);
  if isempty(isspAy); isspAy = ones(size(blk,1),1); end 
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if strcmp(pblk{1},'s')
        n = sum(pblk{2});
        m1 = size(At{p,1},2); 
        if (~isempty(At{p,1}))
           if (ismtpermA)
              tmp = At{p,1}*y(1:m1); 
           else
              tmp = At{p,1}*y(permA(p,1:m1),1); 
           end
           Q{p} = smat(pblk,tmp,isspAy(p));
        else
           Q{p} = sparse(n,n);
        end
        if (length(pblk) > 2) %% for low rank constraints
           len = sum(pblk{3}); 
           m2 = length(pblk{3});
           y2 = y(m1+[1:m2]);
           dd = At{p,3};
           idxD = [0; find(diff(dd(:,1))); size(dd,1)];
           yy2 = mexexpand(diff(idxD),y2); 
           DD = spconvert([dd(:,2:3),dd(:,4).*yy2; len,len,0]);    
           Q{p} = Q{p} + At{p,2}*DD*At{p,2}';
        end
     else
        Q{p} = At{p,1}*y; 
     end 
  end
%%********************************************************* 

