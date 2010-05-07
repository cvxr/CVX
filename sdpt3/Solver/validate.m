%%***********************************************************************
%% validate: validate data
%%
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%***********************************************************************

   function  [blk,At,C,b,dim,nnblk,parbarrier] = ...
              validate(blk,At,C,b,par,parbarrier);

   if (nargin >= 5)
      spdensity = par.spdensity; 
   else
      spdensity = 0.4; 
   end
%%
   if ~iscell(blk); 
      error('validate: blk must be a cell array'); end; 
   if (size(blk,2) < 2)
      error('validate: blk must be a cell array with at least 2 columns');
   end 
   if ~iscell(At) | ~iscell(C); 
      error('validate: At, C must be cell arrays'); 
   end
   if (size(At,1) ~= size(blk,1)) 
      if (size(At,2) == size(blk,1)); 
         At = At'; 
      else      
         error('validate: size of At is not compatible with blk'); 
      end
   end   
   if (size(C,1) ~= size(blk,1)) 
      if (size(C,2) == size(blk,1)) 
         C = C'; 
      else
         error('validate: size of C is not compatible with blk'); 
      end
   end   
   if (min(size(b)) > 1); error('validate: b must be a vector'); end
   if (size(b,1) < size(b,2)); b = b'; end
   m = length(b);
%%
%%-----------------------------------------
%% validate blk, At, C
%%-----------------------------------------
%%
   for p = 1:size(blk,1)
      if (size(blk{p,2},1) > size(blk{p,2},2)) 
         blk{p,2} = blk{p,2}';          
      end
      pblk = blk(p,:); 
      n = sum(pblk{2}); 
      numblk = length(pblk{2});
      if strcmp(pblk{1},'s');
         m1 = size(At{p,1},2); 
         n2 = sum(pblk{2}.*pblk{2});  n22 = sum(pblk{2}.*(pblk{2}+1))/2; 
         ntotal(p) = n22; 
         if ~all(size(C{p}) == n) 
            error('validate: blk and C are not compatible'); 
         end 
         if (norm(C{p}-C{p}',inf) > 1e-13*norm(C{p},inf)); 
            error('validate: C is not symmetric'); 
         end
         if (size(At{p,1}) == [m1, n22] & m1~=n22); 
            At{p,1} = At{p,1}'; 
         end
         if (~isempty(At{p,1})) & (size(At{p,1},1) ~= n22)
            error('validate: blk and At not compatible'); 
         end 
         if (nnz(At{p,1}) < spdensity*n22*m1)
            if ~issparse(At{p,1}); At{p,1} = sparse(At{p,1}); end 
         end
         if (length(pblk) > 2) %% for low rank constraints
            len = sum(pblk{3});
            if (size(pblk{1,3},2) < size(pblk{1,3},1))
                blk{p,3} = blk{p,3}'; 
            end
            if any(size(At{p,2})- [n,len])
               error(' low rank structure specified in blk and At not compatible') 
            end
            if (length(At(p,:)) > 2)  & ~isempty(At{p,3})
               if all(size(At{p,3},2)-[1,4])
                  error(' low rank structure in At{p,3} not specified correctly')
               end
               if (size(At{p,3},2) == 1)
                  if (size(At{p,3},1) < size(At{p,3},2)); At{p,3} = At{p,3}'; end
                  lenn = length(At{p,3});
                  constrnum = mexexpand(pblk{3},[1:length(pblk{3})]');
                  At{p,3} = [constrnum, [1:lenn]', [1:lenn]', At{p,3}];
               elseif (size(At{p,3},2) == 4)
                  dd = At{p,3};
                  [dummy,idxsort] = sort(dd(:,1)); 
                  dd = dd(idxsort,:); 
                  lenn = size(dd,1);
                  idxD = [0; find(diff(dd(:,1))); lenn];
                  ii = zeros(lenn,1); jj = zeros(lenn,1);
                  ss = [0,cumsum(pblk{3})];
                  for k = 1:length(pblk{3})
                     idx = [idxD(k)+1 : idxD(k+1)];
                     ii(idx) = dd(idx,2)+ss(k); %% convert to cumulative indexing
                     jj(idx) = dd(idx,3)+ss(k);
                  end
                  At{p,3} = [dd(:,1),ii,jj,dd(:,4)];
               end
            else
               constrnum = mexexpand(pblk{3},[1:length(pblk{3})]');
               At{p,3} = [constrnum,  [1:len]', [1:len]', ones(len,1)];   
            end
         end
         if (nnz(C{p}) < spdensity*n2) | (numblk > 1); 
            if ~issparse(C{p}); C{p} = sparse(C{p}); end;
         else
            if issparse(C{p}); C{p} = full(C{p}); end; 
         end
      elseif strcmp(pblk{1},'q') | strcmp(pblk{1},'l') | strcmp(pblk{1},'u'); 
         ntotal(p) = n; 
         if (min(size(C{p})) ~= 1 | max(size(C{p})) ~= n); 
            error(['validate: blk and C are not compatible']); 
         end; 
         if (size(C{p},1) < size(C{p},2)); C{p} = C{p}'; end
         if (size(At{p,1}) == [m, n] & m~=n); 
            At{p,1} = At{p,1}'; 
         end
         if ~all(size(At{p,1}) == [n,m]);
            error('validate: blk and At not compatible'); 
         end
         if ~issparse(At{p,1}); 
            At{p,1} = sparse(At{p,1}); 
         end 
         if (nnz(C{p}) < spdensity*n); 
            if ~issparse(C{p}); C{p} = sparse(C{p}); end; 
         else
            if issparse(C{p}); C{p} = full(C{p}); end;
         end;
      else
         error(' blk: some fields are not specified correctly'); 
      end
   end
   if (sum(ntotal) < m) 
      error(' total dimension of C should be > length(b)'); 
   end
%%
%%-----------------------------------------
%% problem dimension
%%-----------------------------------------
%%
   dim = zeros(1,4);  nnblk = zeros(1,2);  
   for p = 1:size(blk,1)
      pblk = blk(p,:);
      if strcmp(pblk{1},'s')
         dim(1) = dim(1) + sum(pblk{2}); 
         nnblk(1) = nnblk(1) + length(pblk{2});
         nn(p) = sum(pblk{2}); 
      elseif strcmp(pblk{1},'q')
         dim(2) = dim(2) + sum(pblk{2}); 
         nnblk(2) = nnblk(2) + length(pblk{2});
         nn(p) = length(pblk{2}); 
      elseif strcmp(pblk{1},'l')
         dim(3) = dim(3) + sum(pblk{2}); 
         nn(p) = sum(pblk{2}); 
      elseif strcmp(pblk{1},'u')
         dim(4) = dim(4) + sum(pblk{2}); 
         nn(p) = sum(pblk{2}); 
      end
   end
%%
%%-----------------------------------------
%% validate parbarrier
%%-----------------------------------------
%% 
    if (nargin == 6)
       if ~iscell(parbarrier); 
	  if (length(parbarrier) == size(blk,1))
	     tmp = parbarrier; 
             clear parbarrier;
	     parbarrier = cell(size(blk,1),1);
	     for p = 1:size(blk,1)
		parbarrier{p} = tmp(p);    
             end
          end
       end 
       if (size(parbarrier,2) > size(parbarrier,1))
          parbarrier = parbarrier'; 
       end
       for p = 1:size(blk,1)
          pblk = blk(p,:); 
          if (size(parbarrier{p},1) > size(parbarrier{p},2))
             parbarrier{p} = parbarrier{p}'; 
          end
          len = length(parbarrier{p}); 
          if strcmp(pblk{1},'s') | strcmp(pblk{1},'q')
             if (len == 1)
                parbarrier{p} = parbarrier{p}*ones(1,length(pblk{2})); 
	     elseif (len == 0)
                parbarrier{p} = zeros(1,length(pblk{2})); 
	     elseif (len ~= length(pblk{2}))
                error('blk and parbarrier not compatible');
             end
          elseif strcmp(pblk{1},'l')
             if (len == 1)
                parbarrier{p} = parbarrier{p}*ones(1,sum(pblk{2})); 
	     elseif (len == 0)
		parbarrier{p} = zeros(1,sum(pblk{2})); 
             elseif (len ~= sum(pblk{2}))
                error('blk and parbarrier not compatible');            
             end
          elseif strcmp(pblk{1},'u')
             parbarrier{p}= zeros(1,sum(pblk{2})); 
          end 
       end
    end
%%***********************************************************************
