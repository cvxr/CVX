%%**********************************************************
%% SDPT3soln_SEDUMIsoln: convert SQLP solution in SDPT3 format to 
%%                       SeDuMi format
%%
%% [xx,yy,zz] = SDPT3soln_SEDUMIsoln(blk,X,y,Z,perm); 
%%
%% usage: load SEDUMI_data_file (containing say, A,b,c,K)
%%        [blk,At,C,b,perm] = read_sedumi(A,b,c,K); 
%%        [obj,X,y,Z]  = sdpt3(blk,At,C,b); 
%%        [xx,yy,zz]   = SDPT3soln_SEDUMIsoln(blk,X,y,Z,perm); 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%**********************************************************

  function [xx,yy,zz] = SDPT3soln_SEDUMIsoln(blk,X,y,Z,perm); 

  yy = y; 
  xx = []; zz = []; 
%%
%% extract unrestricted blk
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if strcmp(pblk{1},'u')
        xx  = [xx; X{p,1}];
        zz  = [zz; Z{p,1}];
     end
  end
%%
%% extract linear blk
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     if strcmp(pblk{1},'l')
        xx  = [xx; X{p,1}];
        zz  = [zz; Z{p,1}];
     end
  end
%%
%% extract second order cone blk 
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if strcmp(pblk{1},'q')  
        xx  = [xx; X{p,1}];
        zz  = [zz; Z{p,1}];
     end
  end
%%
%% extract rotated cone blk 
%%
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if strcmp(pblk{1},'r')  
        xx  = [xx; X{p,1}];
        zz  = [zz; Z{p,1}];
     end
  end
%%
%% extract semidefinite cone blk
%%
  per = []; 
  len = 0; 
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if strcmp(pblk{1},'s')  
        sblk(p) = length(pblk{2}); 
        per = [per, perm{p}]; 
        len = len + sum(pblk{2}.*pblk{2}); 
     end
  end
  sblk = sum(sblk); 
  cnt = 1; 
  Xsblk = cell(sblk,1); Zsblk = cell(sblk,1); 
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if strcmp(pblk{1},'s')
        ss = [0,cumsum(pblk{2})]; 
        numblk = length(pblk{2}); 
	Xp = X{p,1}; 
        Zp = Z{p,1}; 
        xtmp = []; ztmp = []; 
        for tt = 1:numblk
   	   if (numblk > 1) 
              idx = [ss(tt)+1: ss(tt+1)];
              Xsblk{cnt} = full(Xp(idx,idx)); 
              Zsblk{cnt} = full(Zp(idx,idx)); 
           else
              Xsblk{cnt} = Xp; 
              Zsblk{cnt} = Zp; 
           end
	   cnt = cnt + 1; 
        end        
     end
  end
  if ~isempty(per)
     Xsblk(per) = Xsblk; Zsblk(per) = Zsblk; 
     xtmp = zeros(len,1); ztmp = zeros(len,1);
     cnt = 0; 
     for p = 1:sblk
        if strcmp(pblk{1},'s') 
           idx = [1:length(Xsblk{p})^2]; 
           xtmp(cnt+idx) = Xsblk{p}(:); 
           ztmp(cnt+idx) = Zsblk{p}(:); 
           cnt = cnt + length(idx); 
        end
     end
     xx = [xx; xtmp]; 
     zz = [zz; ztmp]; 
  end
%%**********************************************************
