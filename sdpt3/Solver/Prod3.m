%%************************************************************
%% Prod3: compute the entries of Q = A*B*C specified in 
%%        nzlistQ. 
%% 
%% Q = Prod3(blk,A,B,C,sym,nzlistQ)
%% Important: (a) A is assumed to be symmetric if nzlistQ
%%                has 2 columns (since mexProd2nz computes A'*B). 
%%            (b) The 2nd column of nzlistQ must be sorted in 
%%                ascending order. 
%%
%% (optional) sym = 1, if Q is symmetric.
%%                = 0, otherwise. 
%% (optional) nzlistQ = list of non-zero elements of Q to be 
%%                      computed.
%%
%% SDPT3: version 3.1
%% Copyright (c) 1997 by
%% K.C. Toh, M.J. Todd, R.H. Tutuncu
%% Last Modified: 16 Sep 2004
%%************************************************************

  function Q = Prod3(blk,A,B,C,sym,nzlistQ)

  if (nargin<5); sym = 0; end;
  checkcell = [iscell(A) iscell(B) iscell(C)]; 
  if (nargin==6)
     checkcell(1,4) = iscell(nzlistQ); 
  else      
     nzlistQ = inf; 
  end
%%
  if any(checkcell-1)      
     if (size(blk,1) > 1) 
        error('Prod3: blk and A,B,C are not compatible'); 
     end
     if strcmp(blk{1},'s')
        [len,len2] = size(nzlistQ);  
        if (len == 0); nzlistQ = inf; len2 = 1; end; 
        if (len2 == 1) & (nzlistQ == inf) 
           tmp = Prod2(blk,A,B,0);   
           Q = Prod2(blk,tmp,C,sym); 
        else
           tmp = Prod2(blk,B,C,0);
           Q = mexProd2nz(blk,A,tmp,nzlistQ); 
           if sym; Q = 0.5*(Q+Q'); end; 
        end     
     elseif strcmp(blk{1},'q') | strcmp(blk{1},'l') | strcmp(blk{1},'u')
        Q = A.*B.*C;
     end
  else 
     error('Prod3: A,B,C,nzlistQ must all be matrices'); 
  end
%%************************************************************
