%%***************************************************************
%% convertRcone: convert rotated cone to socp cone
%%
%% [blk,At,C,b,T] = convertRcone(blk,At,C,b);
%%
%%***************************************************************
 
  function [blk,At,C,b,T] = convertRcone(blk,At,C,b);

  T = cell(size(blk,1),1); 
  for p = 1:size(blk,1)
     pblk = blk(p,:);      
     if strcmp(pblk{1},'r')
        if (min(pblk{2}) < 3) 
           error('rotated cones must be at least 3-dimensional'); 
        end
        n = sum(pblk{2}); len = length(pblk{2}); 
        ss = [0,cumsum(pblk{2})];
        idx = 1+ss(1:len)'; 
        ir2 = 1/sqrt(2)*ones(len,1);
        dd = [idx,idx,ir2-1; idx,idx+1,ir2; 
              idx+1,idx,ir2; idx+1,idx+1,-ir2-1];
        T{p} = speye(n,n) + spconvert([dd; n,n,0]);   
        blk{p,1} = 'q'; 
        At{p,1} = T{p}*At{p,1}; 
        C{p,1} = T{p}*C{p,1}; 
     end
  end
%%***************************************************************
