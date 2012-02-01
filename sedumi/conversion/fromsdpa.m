function [At,b,c,K]=fromsdpa(fname)
% readsdpa Reads SDP problem from sparse SDPA-formatted input file.
%    [At,b,c,K] = readsdpa(fname) produces conic optimization problem
%    data, as can be used by sedumi. "fname" contains a full pathname
%    to the SDPA 4.10-formatted file. (SDPA is a stand-alone solver for
%    semidefinite programming, by Fujisawa, Kojima and Nakata.  It is
%    used as a standard format in the collection SDPLIB by Borchers.)
%
%    To read and solve the problem "arch0", you may type
%
%    [At,b,c,K] = readsdpa('arch0.dat-s');
%    [x,y,info] = sedumi(At,b,c,K);
%
%    The above 2 lines assume that arch0.dat-s is somewhere in your MATLAB
%    search path, that it is not compressed, and that you know the extension
%    'dat-s'.  To alleviate these conditions, you may like to use the script
%    GETPROBLEM.
%
% SEE ALSO SeDuMi, getproblem, frompack, prelp.

%
% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA%

fid = fopen(fname,'r');
if fid == -1
    error('File not found.')
end
%The number of equality constraints
m='';
while(isempty(m))
    m = fscanf(fid,'%d',1);
    if fgetl(fid)==-1;
        error('Invalid SDPA file.')
    end
end
%The number of blocks
nblocks = fscanf(fid,'%d',1);
if isempty(nblocks)
    error('Invalid SDPA file.')
end
fgetl(fid);
%The dimension of the blocks
%Negative means diagonal block, we convert that to linear
%.,(){} are omitted
dims=sscanf(regexprep(fgetl(fid),'[\.\,(){}]',' '),'%d',nblocks)';
if ~isempty(dims(dims==0)) || length(dims)~=nblocks
    error('Invalid SDPA file.')
end
N=-sum(dims(dims<0))+sum(dims(dims>0).^2);
nblocks=length(dims);
%Create a vector with the offsets of the blocks.
%This is one less than the blockstart!
%Diagonal and one dimensional blocks are coming first
loffset=0;
sdpoffset=sum(abs(dims(dims<=1)));
offset=zeros(1,nblocks);
for i=1:nblocks
    if dims(i)<=1
        offset(i)=loffset;
        loffset=loffset+abs(dims(i));
    else
        offset(i)=sdpoffset;
        sdpoffset=sdpoffset+dims(i)^2;
    end
end
%This is needed so that we can compute where the second column starts
stride=dims;
stride(stride<0)=0;

%Vector b
%,(){} are omitted
b = sscanf(regexprep(fgetl(fid),'[\,(){}]',' '),'%f',m);
if length(b)~=m
    error('Invalid SDPA file.')
end
if nnz(b)/m<0.1
    b=sparse(b);
end

%Coefficients
E = fscanf(fid,'%d %d %d %d %f',[5. inf]);
fclose(fid);

%Extract the objective
cE=E(:,E(1,:)==0);
%Repeating indices in the sparse matrix structure create the sum of
%elements, we need to clear one for the diagonals
data2=cE(5,:);
data2(cE(3,:)==cE(4,:))=0;
c=-sparse([offset(cE(2,:))+(cE(3,:)-1).*stride(cE(2,:))+cE(4,:),offset(cE(2,:))+(cE(4,:)-1).*stride(cE(2,:))+cE(3,:)],1,[cE(5,:),data2],N,1);
clear cE data2
AtE=E(:,E(1,:)~=0);
clear E
data2=AtE(5,:);
data2(AtE(3,:)==AtE(4,:))=0;
At=sparse([offset(AtE(2,:))+(AtE(3,:)-1).*stride(AtE(2,:))+AtE(4,:),offset(AtE(2,:))+(AtE(4,:)-1).*stride(AtE(2,:))+AtE(3,:)],[AtE(1,:),AtE(1,:)],[AtE(5,:),data2],N,m);
clear AtE data2

K.l=-sum(dims(dims<0))+sum(dims==1);
K.s=dims(dims>1);

%Finally, let us correct the sparsity
[i,j,v]=find(c);
c=sparse(i,1,v,N,1);
[i,j,v]=find(At);
At=sparse(i,j,v,N,m);
