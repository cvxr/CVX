%%**********************************************************
%% SDPT3data_SEDUMIdata: convert SQLP data in SDPT3 format to 
%%                       SeDuMi format
%%
%% [At,b,c,K] = SDPT3data_SEDUMIdata(blk,AAt,CC,bb); 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**********************************************************

  function [At,b,c,K] = SDPT3data_SEDUMIdata(blk,AAt,CC,bb); 

  c = []; At = []; 
  b = bb;
  mm = length(bb); 
%%
  if (~iscell(CC))
     Ctmp = CC; clear CC; CC{1} = Ctmp; 
  end
%%
%% extract unrestricted blk
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if (p==1); K.f = []; end      
     if strcmp(pblk{1},'u')
        K.f = [K.f, pblk{2}];
        At  = [At; AAt{p}];
        c   = [c; CC{p}];
     end
  end
  K.f = sum(K.f); 
%%
%% extract linear blk
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if (p==1); K.l = []; end 
     if strcmp(pblk{1},'l')
        K.l = [K.l, pblk{2}];
        At  = [At; AAt{p,1}];
        c   = [c; CC{p,1}];
     end
  end
  K.l = sum(K.l);
%%
%% extract second order cone blk 
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if (p==1); K.q = []; end
     if strcmp(pblk{1},'q')  
        K.q = [K.q, pblk{2}];
        At  = [At; AAt{p,1}];
        c   = [c; CC{p,1}];
     end
  end
%%
%% extract rotated cone blk 
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if (p==1); K.r = []; end
     if strcmp(pblk{1},'r')  
        K.r = [K.r, pblk{2}];
        At  = [At; AAt{p,1}];
        c   = [c; CC{p,1}];
     end
  end
%%
%% extract semidefinite cone blk
%%
  for p = 1:size(blk,1)
     if (p==1); K.s = []; end
     pblk = blk(p,:); 
     if strcmp(pblk{1},'s')
        K.s = [K.s, pblk{2}];
        ss = [0,cumsum(pblk{2})]; 
        idxstart = [0,cumsum(pblk{2}.*pblk{2})]; 
        numblk = length(pblk{2}); 
        nnzA = nnz(AAt{p,1}); 
        II = zeros(2*nnzA,1); 
        JJ = zeros(2*nnzA,1);
        VV = zeros(2*nnzA,1); 
        m2 = size(AAt{p,1},2); 
        if (length(pblk) > 2)
           rr = [0, cumsum(pblk{3})]; 
	   dd = AAt{p,3}; 
           idxD = [0; find(diff(dd(:,1))); size(dd,1)];
        end
        count = 0;      
        for k = 1:mm
   	   if (k<= m2); 
              Ak = smat(pblk,AAt{p,1}(:,k),1); 
           else
              idx = [rr(k)+1 : rr(k+1)];
              Vk  = AAt{p,2}(:,idx);
              len = pblk{3}(k);
              if (size(dd,2) == 4)
                 idx2 = [idxD(k)+1:idxD(k+1)];
                 Dk = spconvert([dd(idx2,2:4); len,len,0]);
              elseif (size(dd,2) == 1); 
 	         Dk = spdiags(dd(idx),0,len,len);  
              end
              Ak = Vk*Dk*Vk';              
           end
           for tt = 1:numblk
              if (numblk > 1)
                 idx = [ss(tt)+1: ss(tt+1)];
                 Aksub = full(Ak(idx,idx)); 
              else
                 Aksub = Ak; 
              end
              tmp = Aksub(:); 
              nzidx = find(tmp); 
              len = length(nzidx);               
              II(count+[1:len],1) = idxstart(tt)+nzidx;
              JJ(count+[1:len],1) = k*ones(length(nzidx),1);
              VV(count+[1:len],1) = tmp(nzidx); 
              count = count + len; 
           end
        end
        II = II(1:count); 
        JJ = JJ(1:count); 
        VV = VV(1:count); 
        At = [At; spconvert([II,JJ,VV; sum(pblk{2}.*pblk{2}), mm, 0])];
	    Cp = CC{p}; 
        ctmp = [];
        for tt = 1:numblk
   	   if (numblk > 1) 
              idx = [ss(tt)+1: ss(tt+1)];
              Csub = full(Cp(idx,idx)); 
           else
              Csub = Cp; 
           end
           ctmp = [ctmp; Csub(:)];
        end
        c = [c; ctmp];
     end
  end
%%**********************************************************
