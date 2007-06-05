%%*******************************************************************
%%  Read in a problem in SDPA sparse format.
%%
%%  [blk,At,C,b] = read_sdpa(fname)
%%
%%  Input: fname = name of the file containing SDP data in
%%                 SDPA foramt. 
%%  Important: the data is assumed to contain only 
%%             semidefinite and linear blocks. 
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%******************************************************************

   function [blk,At,C,b] = read_sdpa(fname); 

%%
%%  Open the file for input
%%
   compressed = 0; 
   if exist(fname)
      fid = fopen(fname,'r');
   elseif exist([fname,'.Z']); 
      compressed = 1; 
      unix(['uncompress ',fname,'.Z']);
      fid = fopen(fname,'r');
   elseif exist([fname,'.gz']); 
      compressed = 2;
      unix(['gunzip ',fname,'.gz']);
      fid = fopen(fname,'r');
   else
      fprintf('*** Problem not found, please specify the correct path or problem name. \n');
      blk = []; At = []; C = []; b = [];
      return;
   end
%%
%%  Clean up special characters and comments from the file 
%%
   [datavec,count] = fscanf(fid,'%c');
   linefeeds = findstr(datavec,char(10));
   comment_chars = '*"=';
   cumidx = [];
   for i=1:length(comment_chars)
      idx = findstr(datavec,comment_chars(i));
      cumidx = [cumidx,idx];
   end
   for j=length(cumidx):-1:1
      if (cumidx(j)==1) | (strcmp(datavec(cumidx(j)-1),char(10)))
         datavec(cumidx(j):linefeeds(min(find(cumidx(j)<linefeeds))))='';
      else
         datavec(cumidx(j):linefeeds(min(find(cumidx(j)<linefeeds)))-1)='';
      end
   end
   special_chars = ',{}()';
   cumidx=[];
   for i=1:length(special_chars)
      idx = findstr(datavec,special_chars(i));
      cumidx = [cumidx,idx];
   end
   datavec(cumidx) = blanks(length(cumidx));
   clear linefeeds;
%%
%% Close the file 
%%
   fclose('all');
   if compressed==1; unix(['compress ',fname]); end;
   if compressed==2; unix(['gzip ',fname]); end; 
%% 
%%  Next, read in basic problem size parameters.
%%
   datavec = sscanf(datavec,'%f'); 
   if size(datavec,1) < size(datavec,2); datavec = datavec'; end; 
   m = datavec(1); 
   numblk  = datavec(2);
   blksize = datavec(2+[1:numblk]); 
   if size(blksize,1) > size(blksize,2); blksize = blksize'; end
%%
%% Get input  b.
%%
   idxstrt = 2+numblk; 
   b = datavec(idxstrt+[1:m]);    
   idxstrt = idxstrt+m; 
   b = -b;
%%
%% Construct blk
%%
   deblksize = 100; 
   spblkidxtmp = find( (blksize>1) & (blksize < deblksize) ); 
   spblkidxtmp = sort(spblkidxtmp);
   deblkidx = find( (blksize<=1) | (blksize >= deblksize) ); 
   denumblk = length(deblkidx); 
   linblkidx = zeros(1,denumblk);  
   for p = 1:denumblk
      n = blksize(deblkidx(p)); 
      if (n > 1); 
         blk{p,1} = 's'; blk{p,2} = n;
         n2 = n*(n+1)/2; 
         At{p,1} = sparse(n2,m);
         C{p,1} = sparse(n,n); 
      else
         linblkidx(p) = p;    
         blk{p,1} = 'l'; blk{p,2} = abs(n); 
         At{p,1} = sparse(abs(n),m); 
         C{p,1} = sparse(abs(n),1); 
      end  
   end
   if ~isempty(spblkidxtmp) 
      maxnumblk = 200; 
      spnumblk = ceil(length(spblkidxtmp)/maxnumblk);
      for q = 1:spnumblk
         if (q < spnumblk)          
            spblkidxall{q} = spblkidxtmp([(q-1)*maxnumblk+1: q*maxnumblk]); 
         else
            spblkidxall{q} = spblkidxtmp([(q-1)*maxnumblk+1: length(spblkidxtmp)]); 
         end
         tmp = blksize(spblkidxall{q}); 
         blk{denumblk+q,1} = 's';  
         blk{denumblk+q,2} = tmp; 
         n2 = sum(tmp.*(tmp+1))/2; 
         At{denumblk+q,1} = sparse(n2,m);
         C{denumblk+q,1} = sparse(sum(tmp),sum(tmp));  
      end
   else
      spnumblk = 0; 
   end
   linblkidx(denumblk+[1:spnumblk]) = zeros(1,spnumblk); 
%%
%% Construct single blocks of A,C
%%
   len = length(datavec);    
   Y = reshape(datavec(idxstrt+1:len),5,(len-idxstrt)/5)';   
   clear datavec;    
   Y = sortrows(Y,[1 2]); 
   matidx = [0; find(diff(Y(:,1)) ~= 0); size(Y,1)];
%%
   for k = 1:length(matidx)-1
      idx = [matidx(k)+1 : matidx(k+1)];       
      Ytmp  = Y(idx,1:5); 
      matno = Ytmp(1,1); 
      Ytmp2 = Ytmp(:,2); 
      for p = 1:denumblk 
         n  = blksize(deblkidx(p));   
         idx = find(Ytmp2 == deblkidx(p)); 
         ii = Ytmp(idx,3); jj = Ytmp(idx,4); vv =Ytmp(idx,5); 
         len = length(idx); 
         if (n > 1)
            idxtmp = find(ii > jj); 
            if ~isempty(idxtmp); 
               tmp = jj(idxtmp); 
               jj(idxtmp) = ii(idxtmp); ii(idxtmp) = tmp; 
            end
            tmp = -sparse(ii,jj,vv,n,n); 
            tmp = tmp + triu(tmp,1)'; 
         else
            tmp = -sparse(ii,ones(len,1),vv,abs(n),1); 
         end
         if (matno == 0) 
            C{p,1} = tmp; 
         else
            if (n > 1)
               At{p,1}(:,matno) = svec(blk(p,:),tmp,1); 
            else
               At{p,1}(:,matno) = tmp; 
            end
         end
      end
   end 
%%
%% Construct big sparse block of A,C 
%%
if (spnumblk > 0)
   Y1 = Y(:,1); 
   diffY1 = find(diff([-1; Y1; inf]));  
   for kk = 1:length(diffY1)-1
      idx = [diffY1(kk) : diffY1(kk+1)-1];
      matno = Y1(diffY1(kk)); 
      Ytmp  = Y(idx,1:5);
      Ytmp2 = Ytmp(:,2); 
      maxYtmp2  = Ytmp2(length(Ytmp2)); 
      minYtmp2  = Ytmp2(1);
      diffYtmp2 = Ytmp2(find(diff([-1; Ytmp2]))); 
      for q = 1:spnumblk
         spblkidx    = spblkidxall{q};
         maxspblkidx = spblkidx(length(spblkidx));
         minspblkidx = spblkidx(1); 
         count = 0; 
         if (minYtmp2 <= maxspblkidx) & (maxYtmp2 >= minspblkidx)
            tmpblksize = blksize(spblkidx);
            n = sum(tmpblksize); 
            cumspblksize = [0 cumsum(tmpblksize)]; 
            n2 = sum(tmpblksize.*(tmpblksize+1))/2; 
            idx = zeros(n2,1); offset = zeros(n2,1); 
            for t = [1:length(diffYtmp2)] 
  	       p = find(spblkidx == diffYtmp2(t)); 
               if ~isempty(p)
	          idxtmp = find(Ytmp2 == spblkidx(p));
                  len = length(idxtmp); 
     	          idx(count+[1:len]) = idxtmp; 
                  offset(count+[1:len]) = cumspblksize(p); 
                  count = count + len; 
               end
            end 
            idx = idx(1:count); offset = offset(1:count); 
            ii = Ytmp(idx,3)+offset; jj = Ytmp(idx,4)+offset; vv = Ytmp(idx,5);
            idxtmp = find(ii > jj); 
            if ~isempty(idxtmp); 
               tmp = jj(idxtmp); 
               jj(idxtmp) = ii(idxtmp); ii(idxtmp) = tmp; 
            end
	    idxeq = find(ii==jj); 
            tmp = spconvert([ii jj -vv; jj ii -vv; n n 0]) ...
                  + spconvert([ii(idxeq) jj(idxeq) vv(idxeq); n n 0]);
            if (matno == 0) 
               C{denumblk+q,1} = tmp;
            else
               At{denumblk+q,1}(:,matno) = svec(blk(denumblk+q,:),tmp,1); 
            end
         end
      end
   end
end
%%
%% put all linear blocks together as a single linear block
%% 
   idx = find(linblkidx); 
   if (length(idx) > 1)
      sdpidx = find(linblkidx==0); 
      blktmp = 0; Atmp = []; Ctmp = [];
      for k = 1:length(idx)
          tmp = linblkidx(idx(k)); 
          blktmp = blktmp+blk{tmp,2}; 
          Atmp = [Atmp; At{tmp}];
          Ctmp = [Ctmp; C{tmp}]; 
      end
      At = At(sdpidx); C = C(sdpidx); blk = blk(sdpidx,:);    
      len = length(sdpidx); 
      blk(2:len+1,:) = blk; 
      blk{1,1} = 'l'; blk{1,2} = blktmp; 
      At(2:len+1,1) = At; C(2:len+1,1) = C; 
      At{1,1} = Atmp; C{1,1} = Ctmp;   
   end
%%******************************************************************

