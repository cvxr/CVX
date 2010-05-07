%%********************************************************************
%% infeaspt: generate an initial point for sdp.m
%%
%%  [X0,y0,Z0] = infeaspt(blk,At,C,b,options,scalefac);
%%
%%  options = 1  if want X0,Z0 to be scaled identity matrices
%%          = 2  if want X0,Z0 to be scalefac*(identity matrices).
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%********************************************************************

   function [X0,y0,Z0] = infeaspt(blk,At,C,b,options,scalefac);
%%
   if (nargin < 5); options = 1; end;
   if (options == 1); scalefac = []; end;
   if (options == 2) & (nargin < 6); scalefac = 1000; end;
   if (scalefac <= 0); error('scalefac must a positive number'); end;
   state = rand('state'); 
   rand('state',0);
%%
   if ~iscell(At); At = {At}; end;
   if ~iscell(C);  C = {C}; end;
   m = length(b);       
   if all(size(At) == [size(blk,1) m]); 
      convertyes = zeros(size(blk,1),1); 
      for p = 1:size(blk,1)
         if strcmp(blk{p,1},'s') & all(size(At{p,1}) == sum(blk{p,2}))
            convertyes(p) = 1;    
         end
      end
      if any(convertyes)
         At = svec(blk,At,ones(size(blk,1),1));
      end
   end; 
%%
   %%[blk,At,C,b] = validate(blk,At,C,b);
%%
   X0 = cell(size(C)); Z0 = cell(size(C));
   m = length(b); 
   for p = 1:size(blk,1); 
      pblk = blk(p,:); 
      blktmp = pblk{2};
      n = length(C{p});
      y0 = zeros(m,1);
      b2 = 1 + abs(b');
      if (options == 1);
         if strcmp(pblk{1},'s');
            normAni = [];
            X0{p} = sparse(n,n); Z0{p} = sparse(n,n);
            ss = [0, cumsum(blktmp)];
            tt = [0, cumsum(blktmp.*(blktmp+1)/2)];
            for i = 1:length(pblk{2})
               if ~isempty(At{p,1})
                  pos = [tt(i)+1 : tt(i+1)];
                  Ai = At{p,1}(pos,:);
                  normAni = 1+sqrt(sum(Ai.*Ai));
               end
               if (length(At(p,:)) >= 2) %% for low rank constraints
                  dd = At{p,3};
                  qq = [0, cumsum(pblk{3})]; normtmp = ones(1,length(pblk{3}));
                  idxD = [0; find(diff(dd(:,1))); size(dd,1)];
                  for k = 1:length(pblk{3})
                      idx = [qq(k)+1 : qq(k+1)];
                      idx2 = [idxD(k)+1: idxD(k+1)];
                      Ak = At{p,2}(:,idx);
                      ii = dd(idx2,2)-qq(k); %% undo cumulative indexing 
                      jj = dd(idx2,3)-qq(k);
                      len = pblk{3}(k);
                      Dk = spconvert([ii,jj,dd(idx2,4); len,len,0]);
                      tmp = Ak'*Ak*Dk;
                      normtmp(1,k) = 1+sqrt(sum(sum(tmp.*tmp'))); 
                  end
                  normAni = [normAni, normtmp];
               end
               pos = [ss(i)+1 : ss(i+1)];  ni = length(pos);
               tmp = C{p}(pos,pos);
               normCni = 1+sqrt(sum(sum(tmp.*tmp)));
               const  = 10; %%--- old: const = 1; 
               constX = max([const,sqrt(ni),ni*(b2./normAni)]); 
               constZ = max([const,sqrt(ni),normAni,normCni]);
               X0{p}(pos,pos) = constX*spdiags(1+1e-10*rand(ni,1),0,ni,ni);
               Z0{p}(pos,pos) = constZ*spdiags(1+1e-10*rand(ni,1),0,ni,ni);
            end
         elseif strcmp(pblk{1},'q');
            s = 1+[0, cumsum(blktmp)];
            len = length(blktmp);
            normC = 1+norm(C{p});
            normA = 1+sqrt(sum(At{p,1}.*At{p,1}));
            idenqX = zeros(sum(blktmp),1);
            idenqZ = zeros(sum(blktmp),1);
            idenqX(s(1:len)) = max([1,b2./normA])*sqrt(blktmp') ;
            idenqZ(s(1:len)) = max([sqrt(blktmp); max([normA,normC])*ones(1,len)])';
            idenqX(s(1:len)) = idenqX(s(1:len)).*(1+1e-10*rand(len,1)); 
            idenqZ(s(1:len)) = idenqZ(s(1:len)).*(1+1e-10*rand(len,1)); 
            X0{p} = idenqX;
            Z0{p} = idenqZ;
         elseif strcmp(pblk{1},'l');
            normC = 1+norm(C{p});
            normA = 1+sqrt(sum(At{p,1}.*At{p,1}));
            const = 10; %%--- old: const =1; 
            constX = max([const,sqrt(n),sqrt(n)*b2./normA]); 
            constZ = max([const,sqrt(n),normA,normC]);
            X0{p} = constX*(1+1e-10*rand(n,1));
            Z0{p} = constZ*(1+1e-10*rand(n,1));
         elseif strcmp(pblk{1},'u');
            X0{p} = sparse(n,1);
            Z0{p} = sparse(n,1);
         else
            error(' blk: some fields not specified correctly'); 
         end;
      elseif (options == 2);
         if strcmp(pblk{1},'s');
            n = sum(blktmp); 
            X0{p} = scalefac*spdiags(1+1e-10*rand(n,1),0,n,n); 
            Z0{p} = scalefac*spdiags(1+1e-10*rand(n,1),0,n,n); 
         elseif strcmp(pblk{1},'q');
            s = 1+[0, cumsum(blktmp)];
            len = length(blktmp);
            idenq = zeros(sum(blktmp),1);
            idenq(s(1:len)) = 1+1e-10*rand(len,1);
            X0{p} = scalefac*idenq;
            Z0{p} = scalefac*idenq;
         elseif strcmp(pblk{1},'l');
            X0{p} = scalefac*(1+1e-10*rand(n,1));
            Z0{p} = scalefac*(1+1e-10*rand(n,1));
         elseif strcmp(pblk{1},'u');
            X0{p} = sparse(n,1);
            Z0{p} = sparse(n,1);
         else
            error(' blk: some fields not specified correctly'); 
         end
      end
   end
   rand('state',state); 
%%********************************************************************
