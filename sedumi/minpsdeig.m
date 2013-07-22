%                                                mineig = minpsdeig(x,K)
% MINPSDEIG  Computes the smallest spectral coefficients of x w.r.t. K
% Uses an iterative method if the matrix is large, takes the minimum of all
% the eigenvalues if the matrix is small.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also sedumi

function mineig = minpsdeig(x,K)
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

% disp('The SeDuMi binaries are not installed.')
% disp('In Matlab, launch "install_sedumi" in the folder you put the SeDuMi files.')
% disp('For more information see the file Install.txt.')
% error(' ')

% [labold,qold]=psdeig(x,K);
Ks=K.s;
if isempty(Ks)
    mineig=[];
    return
end

startindices=K.sblkstart-K.mainblks(end)+1;
ncones=length(Ks);
OPTS.disp=0;
mineigvector=zeros(ncones,1);
for k = 1:ncones
    Xk = reshape(x(startindices(k):startindices(k+1)-1),Ks(k),Ks(k));
    if Ks(k)>500
        mineigvector(k) = eigs(Xk + Xk',1,'SA',OPTS);
    else
        mineigvector(k) = min(eig(Xk + Xk'));
    end
end
mineig = min(mineigvector) / 2;
