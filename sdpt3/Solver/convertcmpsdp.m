%%*********************************************************
%% convertcmpsdp: convert SDP with complex data into one
%%                with real data by converting 
%% 
%%  C - sum_{k=1}^m yk*Ak psd
%%  to 
%%  [CR,-CI] - sum ykR*[AkR,-AkI] psd
%%  [CI, CR]           [AkI, AkR]  
%%
%%  ykI = 0 for k = 1:m
%%
%% [bblk,AAt,CC,bb] = convertcmpsdp(blk,A,C,b);
%%
%%*********************************************************

  function [bblk,AAt,CC,bb,iscmp] = convertcmpsdp(blk,A,C,b);

  m = length(b); 
  [pp,mm] = size(A); 
  if (pp ~= size(blk,1)) 
     error('blk and A not compatible'); 
  end
  numblk = size(blk,1); 
  iscmp = zeros(numblk,m+1); 
  for p = 1:size(blk,1)
     pblk = blk(p,:);
     len = size(A(p),2); 
     for k = 1:len
        if ~isempty(A{p,k})
           iscmp(p,k) = 1-isreal(A{p,k}); 
        end
     end
     iscmp(p,m+1) = 1-isreal(C{p});
  end
  iscmp = norm(iscmp,'fro');
%%
  if (iscmp == 0) 
     %% data is real
     bblk = blk; AAt = A; CC = C; bb = b; 
     return; 
  end
%%
  bb = [real(b)]; 
  bblk = cell(size(blk,1),2); 
  for p = 1:size(blk,1)
     pblk = blk(p,:); 
     if (size(pblk{2},1) > size(pblk{2},2))
        pblk{2} = pblk{2}'; 
     end
     if strcmp(pblk{1},'s')
        ss = [0,cumsum(pblk{2})]; 
        ss2 = [0,cumsum(2*pblk{2})];
        n = sum(pblk{2}); 
        n2 = sum(pblk{2}.*(pblk{2}+1))/2; 
        AR = cell(1,m); Ctmp = sparse(2*n,2*n);  
        if (size(A{p},1)==n2 & size(A{p},2)==m); 
           Atype = 1; 
        elseif (size(A(p),1)==1 & size(A(p),2)==1); 
           Atype = 2; 
        else
           error('convertcmp: At is not properly coded');  
        end
        for k = 0:m
           if (k == 0)
              Ak = C{p}; 
           else
              if (Atype == 1)
                 Ak = smat(pblk,A{p}(:,k),1); 
              elseif (Atype == 2) 
                 Ak = A{p,k}; 
              end
           end
           Atmp = sparse(2*n,2*n);
           if (length(pblk{2}) == 1)
              tmp = [real(Ak),-imag(Ak);  imag(Ak), real(Ak)];  
              if (k==0)
                 Ctmp = tmp;
              else
                 Atmp = tmp;
              end
           else            
              for j = 1:length(pblk{2})
                 idx = [ss(j)+1: ss(j+1)];              
                 Akj = Ak(idx,idx); 
                 tmp = [real(Akj),-imag(Akj);  imag(Akj), real(Akj)]; 
                 idx2 = [ss2(j)+1: ss2(j+1)]; 
                 if (k==0)
                    Ctmp(idx2,idx2) = tmp; 
                 else
   	            Atmp(idx2,idx2) = tmp;  
                 end 
              end
           end
           if (k==0);
              CC{p,1} = Ctmp; 
           else
              AR{k} = Atmp; 
           end
        end
        bblk{p,1} = 's'; bblk{p,2} = 2*pblk{2}; 
        AAt(p,1) = svec(bblk(p,:),AR); 
     elseif strcmp(pblk{1},'q'); 
        error('SOCP block with complex data is currently not allowed'); 
     elseif strcmp(pblk{1},'l');
        if isreal(A{p}) & isreal(C{p})
           bblk(p,:) = blk(p,:); 
           AAt{p,1} = A{p}; CC{p,1} = C{p}; 
        else
           error('data for linear block must be real'); 
        end
     elseif strcmp(pblk{1},'u');  
        if isreal(A{p}) & isreal(C{p})
           bblk(p,:) = blk(p,:); 
           AAt{p,1} = A{p}; CC{p,1} = C{p}; 
        else
           error('data for unrestricted block must be real'); 
        end
     end      
  end
%%*********************************************************
