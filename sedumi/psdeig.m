%                                                [lab,q] = psdeig(x,K)
% PSDEIG  Computes spectral coefficients of x w.r.t. K
%   Arguments "q" is optional - without it's considerably faster.
%   FLOPS indication: 1.3 nk^3 versus 9.0 nk^3 for nk=500,
%                     1.5 nk^3        9.8 nk^3 for nk=50.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function [lab,q] = psdeig(x,K) %#ok
%
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
%

Ks=K.s;
if isempty(Ks)
    lab=[];
    return
end

lab(sum(Ks),1)=0;
startindices=K.sblkstart-K.mainblks(end)+1;
labindices=cumsum([1,Ks]);
ncones=length(Ks);
if nargout==1
    %only eigenvalues are needed, not eigenvectors
    for k = 1:ncones
        Xk = reshape(x(startindices(k):startindices(k+1)-1),Ks(k),Ks(k));
        lab(labindices(k):labindices(k+1)-1) = eig(Xk + Xk');
    end
else
    %eigenvalues and eigenvectors
    q=zeros(sum(Ks.^2),1);
    for k = 1:ncones
        Xk = reshape(x(startindices(k):startindices(k+1)-1),Ks(k),Ks(k));
        [q(startindices(k):startindices(k+1)-1), temp] = eig(Xk + Xk');
        lab(labindices(k):labindices(k+1)-1)=diag(temp);
    end
end
%We actually got the eigenvalues of twice the matrices, so we need to scale
%them back.
lab=lab/2;
