function [newAt,newb,newc,newK,prepinfo]=preprocessSDP(At,b,c,K)
%[newAt,newb,newc,newK,prepinfo]=preprocessSDP(At,b,c,K)
%Preprocess the SDP part of a problem, return the new variables and the
%info needed to postprocess the solutions at the end. 
%
%prepinfo: a cell array describing the operations for the cones
%prepinfo{i}(1)=0: No change, prepinfo{i}(2) is the number of columns that are
%                  unchanged.
%              =1: Diagonal matrix converted to linear variables,
%                  prepinfo{i}(2) is the number of linear variables that are
%                  extracted.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi, postprocessSDP, pretransfo

% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
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
% 02110-1301, USA

blockstart=K.l+K.f +sum(K.q)+[0;reshape(cumsum((K.s).^2),length(K.s),1)]+1;
sparsepattern=c-At*sparse(randn(length(b),1));

newAt=[];
newc=[];
newK=K;
newK.s=[];
newb=b;
prepinfo={[0; K.l + K.f + sum(K.q)]};

for i=1:length(K.s)
    %Get aggregate sparsity pattern for the ith block
    blocksparsepattern=reshape(sparsepattern(blockstart(i):blockstart(i+1)-1),K.s(i),K.s(i));
    % Uncomment the following line if you want to see the sparsity pattern
    % of the blocks in s
    %figure,spy(blocksparsepattern)

    %Extract information from the sparsity pattern
    [rindex,cindex]=find(blocksparsepattern);
    clear blocksparsepattern

    %Calculate bandwidth (0 means diagonal matrix)
    bandwidth=norm(rindex-cindex,Inf);
    clear rindex cindex;
    if bandwidth==0
        %We have a diagonal block, we convert it to nonnegative variables.
        newK.l=newK.l+K.s(i);
        prepinfo{end+1}=[1; K.s(i)];
    else
        newK.s=[newK.s;K.s(i)];
        %If the previous block is also unchanged then we move them together.
        if  prepinfo{end}(1)==0
            prepinfo{end}(2)=prepinfo{end}(2)+K.s(i)^2;
        else
            prepinfo{end+1}=[0; K.s(i)^2];
        end
    end
end


%Now we transform At and c. We do this here to gain speed if there are a
%lot of small cones.
for j=1:length(prepinfo)
    op=prepinfo{j};
    switch op(1)
        case 0
            %Do nothing
            newAt=[newAt;At(1:op(2),:)];
            newc=[newc;c(1:op(2))];
            At=At(op(2)+1:end,:);
            c=c(op(2)+1:end);
        case 1
            %Diagonal matrix
            newAt=[[At(1:op(2)+1:op(2)^2,:),sparse(op(2),size(newAt,2)-size(At,2))];newAt];
            newc=[c(1:op(2)+1:op(2)^2);newc];
            At=At(op(2)^2+1:end,:);
            c=c(op(2)^2+1:end);
            newK.rsdpN=newK.rsdpN-1;
    end
end

