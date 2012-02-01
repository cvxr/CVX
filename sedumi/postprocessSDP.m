function [x,y,K]=postprocessSDP(newx,newy,prepinfo,newK)
% [x,y,K]=postprocessSDP(newx,newy,prepinfo,newK)
% postprocessSDP: Postprocesses an SDP solution using the info from
% preprocessing (preprocessSDP)
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi, preprocessSDP

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

if issparse(newx)
    x=sparse([]);
else
    x=[];
end
y=newy;
K=newK;

for j=length(prepinfo):-1:1
    op=prepinfo{j};
    switch op(1)
        case 0
            %Do nothing
            x=[newx(end-op(2)+1:end);x];
            newx=newx(1:end-op(2));
        case 1
            %Diagonal matrix
            x=[reshape(diag(newx(1:op(2))),op(2)^2,1);x];
            newx=newx(op(2)+1:end);
            K.l=K.l-op(2);
            K.s(end+1)=op(2);
            K.rsdpN=K.rsdpN+1;
    end
end
