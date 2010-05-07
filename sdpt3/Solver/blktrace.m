%%**********************************************************************
%% blktrace: compute <X1,Z1> + ... + <Xp,Zp>
%%              
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**********************************************************************

  function trXZ = blktrace(blk,X,Z,parbarrier); 

  if (nargin == 3) 
     trXZ = 0; 
     for p = 1:size(blk,1)
        pblk = blk(p,:);
        if strcmp(pblk{1},'s')
           if (length(pblk{2}) == 1)
              trXZ = trXZ + sum(sum(X{p}.*Z{p})); 
           else
              xx = mexsvec(pblk,X{p},0); 
              zz = mexsvec(pblk,Z{p});
              trXZ = trXZ + xx'*zz;
           end
        else
           trXZ = trXZ + sum(X{p}.*Z{p}); 
        end
     end
  elseif (nargin == 4)
     trXZ = 0; 
     for p = 1:size(blk,1)
        pblk = blk(p,:);
        if (norm(parbarrier{p}) == 0)
           if strcmp(pblk{1},'s')
              if (length(pblk{2}) == 1)
                 trXZ = trXZ + sum(sum(X{p}.*Z{p})); 
              else
                 xx = mexsvec(pblk,X{p},0); 
                 zz = mexsvec(pblk,Z{p});
                 trXZ = trXZ + xx'*zz;
              end
           else
              trXZ = trXZ + sum(X{p}.*Z{p}); 
           end
        else
           idx = find(parbarrier{p} == 0); 
           if ~isempty(idx)             
              if strcmp(pblk{1},'s')
                 sumXZ = sum(X{p}.*Z{p}); 
                 ss = [0,cumsum(pblk{2})];
                 for k = 1:length(idx)
                    idxtmp = [ss(idx(k))+1:ss(idx(k)+1)];
                    trXZ = trXZ + sum(sumXZ(idxtmp)); 
                 end
              elseif strcmp(pblk{1},'q')
	         tmp = qops(pblk,X{p},Z{p},1);
                 trXZ = trXZ + sum(tmp(idx));
              elseif strcmp(pblk{1},'l')
                 trXZ = trXZ + sum(X{p}(idx).*Z{p}(idx)); 
              end
           end
        end
     end
  end
%%**********************************************************************




