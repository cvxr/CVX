%  This function takes a problem in SeDuMi MATLAB format and writes it out 
%  in SDPpack format.  
%
%  Usage:
%
%  writesdp(fname,A,b,c,K)
%
%      fname           Name of SDPpack file, in quotes
%      A,b,c,K         Problem in SeDuMi form
%
%  Notes:
%
%     Problems with complex data are not allowed.    
%
%     Rotated cone constraints are not supported.  
%
%     Nonsymmetric A.s and C.s matrices are symmetrized with A=(A+A')/2
%     a warning is given when this happens.
%
%     Floating point numbers are written out with 18 decimal digits for
%     accuracy.
%
%  Please contact the author (Brian Borchers, borchers@nmt.edu) with any
%  questions or bug reports.
% 
function writesdp(fname,A,b,c,K)

%From: 	Brian Borchers[SMTP:borchers@nmt.edu]
%Sent: 	Wednesday, December 08, 1999 12:33 PM
%To: 	J.Sturm@KE.UNIMAAS.NL
%
%Here's a MATLAB routine that will take a problem in SeDuMi's MATLAB format
%and write it out in SDPpack format.

%
%  First, check for complex numbers in A, b, or c.
%
if (isreal(A) ~= 1),
  disp('A is not real!');
  return;
end;
if (isreal(b) ~= 1),
  disp('b is not real!');
  return;
end;
if (isreal(c) ~= 1),
  disp('c is not real!');
  return;
end;
%
%  Check for any rotated cone constraints.
%
if (isfield(K,'r') && (~isempty(K.r)) && (K.r ~= 0)),
  disp('rotated cone constraints are not yet supported.');
  return;
end; 
%
%  Get the size data.
%
if (isfield(K,'l')),
  nlin=K.l;
  sizelin=nlin;
  if (isempty(sizelin)),
    sizelin=0;
    nlin=0;
  end;
else
  nlin=0;
  sizelin=0;
end;

if (isfield(K,'s')),
  nsdpblocks=length(K.s);
  %sizesdp=sum((K.s).^2);
  %if (isempty(sizesdp)),
    %sizesdp=0;
    %nsdpblocks=0;
  %end;
else
  %sizesdp=0;
  nsdpblocks=0;
end;

if (isfield(K,'q')),
  nqblocks=length(K.q);
  sizeq=sum(K.q);
  if (isempty(sizeq)),
    sizeq=0;
    nqblocks=0;
  end;
else
  nqblocks=0;
  sizeq=0;
end;

m=length(b);

%
%  Open up the file for writing.
%
fid=fopen(fname,'w');
%
%  Print out m, the number of constraints.
%
fprintf(fid,'%d \n',m);
%
%  Next, b, with one entry per line.
%
fprintf(fid,'%.18e\n',full(b));
%
%  Next, the semidefinite part.
%
if (nsdpblocks == 0),
  fprintf(fid,'0\n');
else
%
%  Print out the number of semidefinite blocks.
%
  fprintf(fid,'%d\n',nsdpblocks);
%
%  For each block, print out its size.
%
  fprintf(fid,'%d\n',full(K.s));
%
%  Next, the cost matrix C.s.
%
%
%  First, calculate where in c things start.
%
  base=sizelin+sizeq+1;
%
%  Next, work through the blocks.
%
  for i=1:nsdpblocks,
    fprintf(fid,'1\n');
    work=c(base:base+K.s(i)^2-1);
    work=reshape(work,K.s(i),K.s(i));
    if nnz(work ~= work'),
      disp('Non symmetric C.s matrix!');
      work=(work+work')/2;
    end;
    work=triu(work);
    [II,JJ,V]=find(work);
    cnt=length(II);
    fprintf(fid,'%d\n',cnt);
    if (cnt ~= 0),
      fprintf(fid,'%d\n%d\n%.18e\n',[II JJ V]');
    end;
%
%  Next, update to the next base.
%
    base=base+K.s(i)^2;
  end;
%
%  Now, loop through the constraints, one at a time.
%
  for cn=1:m,
%
%  Print out the SDP part of constraint cn.
%
    base=sizelin+sizeq+1;
    for i=1:nsdpblocks,
      fprintf(fid,'1\n');
      work=A(cn,base:base+K.s(i)^2-1);
      work=reshape(work,K.s(i),K.s(i));
      if nnz(work ~= work'),
        disp('Non symmetric A.s matrix!');
        work=(work+work')/2;
      end;
      work=triu(work);

      [II,JJ,V]=find(work);
      cnt=length(II);
      fprintf(fid,'%d\n',cnt);
      if (cnt ~= 0),
        fprintf(fid,'%d\n%d\n%.18e\n',[II JJ V]');
      end;
%
%  Next, update to the next base.
%
      base=base+K.s(i)^2;
    end;
%
% Done with constraint cn
%
  end;
%
% Done with SDP part.
%
end;
%
% Next, handle the Quadratic part.
%
%
% Describe the Q blocks.
%
if (nqblocks == 0),
  fprintf(fid,'0\n');
else
  fprintf(fid,'%d\n',nqblocks);
  fprintf(fid,'%d\n',full(K.q));
%
%  Find C.q.
%
  base=sizelin+1;
  cq=c(base:base+sizeq-1);
%
% Print out the C.q coefficients.
%
  fprintf(fid,'%.18e\n',full(cq));
%
%  Next, the constraint matrix A.q. 
%
  Aq=A(:,base:base+sizeq-1);
%
%  Print out the count of nonzeros.
%
  [II,JJ,V]=find(Aq);
  cnt=length(II);
  fprintf(fid,'1\n');
  fprintf(fid,'%d\n',cnt);
  if (cnt ~= 0),
    fprintf(fid,'%d\n%d\n%.18e\n',[II JJ V]');
  end;
%
% End of handling quadratic part.
%
end;
%
%
% Finally, handle the linear part.
%
if (nlin == 0),
  fprintf(fid,'0\n');
else
%
% Print out the number of linear variables.
%
  fprintf(fid,'%d\n',nlin);
%
% Print out C.l
%
  fprintf(fid,'%.18e\n',full(c(1:nlin)));
%
%  Print out the A matrix.
%
  Al=A(:,1:nlin);
  [II,JJ,V]=find(Al);
  cnt=length(II);
  fprintf(fid,'1\n');
  fprintf(fid,'%d\n',cnt);
  if (cnt ~= 0),
    fprintf(fid,'%d\n%d\n%.18e\n',[II JJ V]');
  end;
end;






