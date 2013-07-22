%                                                    z = psdmul(x,y, K)
% PSDMUL  for full x,y. Computes (XY+YX)/2
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function z = psdjmul(x,y, K)
%
% This file is part of SeDuMi 1.3 by Imre Polik
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
%

startindices=K.sblkstart;
z=zeros(K.N-startindices(1)+1,1);
Ks=K.s;
endindex=K.mainblks(end);
clear K
for k = 1:length(Ks)
    Ksk=Ks(k);
    startk=startindices(k);
    startkp=startindices(k+1);
    temp=reshape(x(startk:startkp-1),Ksk,Ksk)*reshape(y(startk:startkp-1),Ksk,Ksk);
    z(startk-endindex+1:startkp-endindex)=temp+temp';
end
z=z/2;