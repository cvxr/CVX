%%*******************************************************************
%% combine_blk: combine small SDP blocks together, 
%%              combine all SOCP blocks together, etc
%%
%% [blk2,At2,C2,blkinfo] = combine_blk(blk,At,C); 
%%
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%*******************************************************************

   function [blk2,At2,C2,blkinfo] = combine_blk(blk,At,C); 

   blkinfo = zeros(size(blk,1),1);  
   for p = 1:size(blk,1)
      pblk = blk(p,:); 
      if strcmp(pblk{1},'s')
         if (sum(pblk{2}) < 100)
            blkinfo(p) = 1; 
         end
      elseif strcmp(pblk{1},'q')
         blkinfo(p) = 2; 
      elseif strcmp(pblk{1},'r')
         blkinfo(p) = 3;
      elseif strcmp(pblk{1},'l')
         blkinfo(p) = 4;  
      elseif strcmp(pblk{1},'u')
         blkinfo(p) = 5; 
      end
   end
   numblk0 = length(find(blkinfo == 0)); 
   numblk = numblk0 + length(union(blkinfo(blkinfo > 0),[])); 
   blk2 = cell(numblk,2); At2 = cell(numblk,1); C2 = cell(numblk,1); 
   cnt = 0; 
   idx = find(blkinfo==0); %% larger SDP blocks
   if ~isempty(idx)
      len = length(idx); 
      blk2(1:len,:) = blk(idx,:); 
      At2(1:len) = At(idx); C2(1:len) = C(idx);  
      cnt = len; 
   end
   idx = find(blkinfo==1); %% smaller SDP blocks
   Ctmp = []; idxstart = 0; 
   if ~isempty(idx)
      cnt = cnt + 1; 
      blk2{cnt,1} = 's'; blk2{cnt,2} = [];  
      len = length(idx); 
      for k = 1:len
         blk2{cnt,2} = [blk2{cnt,2}, blk{idx(k),2}]; 
         At2{cnt} = [At2{cnt}; At{idx(k)}];
         [ii,jj,vv] = find(C{idx(k)});  
         Ctmp = [Ctmp; [idxstart+ii,idxstart+jj,vv]]; 
         idxstart = idxstart + sum(blk{idx(k),2}); 
      end
   end
   n = sum(blk2{cnt,2}); 
   C2{cnt} = spconvert([Ctmp; n,n,0]);
%%
   for L = [2:5]
      idx = find(blkinfo==L); 
      if ~isempty(idx)
         cnt = cnt + 1; 
         if (L==2) 
            blk2{cnt,1} = 'q';
	 elseif (L==3)
            blk2{cnt,1} = 'r';
	 elseif (L==4)
            blk2{cnt,1} = 'l';
	 elseif (L==5)
            blk2{cnt,1} = 'u';
         end
         blk2{cnt,2} = [];  
         len = length(idx); 
         for k = 1:len
            blk2{cnt,2} = [blk2{cnt,2}, blk{idx(k),2}]; 
            At2{cnt} = [At2{cnt}; At{idx(k)}];
            C2{cnt} = [C2{cnt}; C{idx(k)}]; 
         end
      end
   end
%%*******************************************************************
