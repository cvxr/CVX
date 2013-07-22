%                                           y = psdinvscale(ud,x,K ,transp)
% PSDINVSCALE  Computes length lenud (=sum(K.s.^2)) vector y.
%    Computes y = D(d^{-1}) x with d in K.
%    Y = Ud' \ X / Ud
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also scaleK, factorK.

function y = psdinvscale(ud,x,K)

% This file is part of SeDuMi 1.3 by Imre Polik and Oleksandr Romanko
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

Ks=K.s;
if ~isempty(Ks)
    N=sum(Ks.^2);
    y(N,1)=0;
    startindices=K.sblkstart-K.mainblks(end)+1;
    %Sometimes x containts only the PSD part, sometimes the whole thing
    xstartindices=startindices+(length(x)-N);
    for i=1:K.rsdpN
        Ksi=Ks(i);
        Ksi2=Ksi^2;
        temp=triu(reshape(ud(startindices(i):startindices(i+1)-1),Ksi,Ksi));
        if nnz(temp)<0.05*Ksi2;
            temp=sparse(temp);
        end
        X=reshape(x(xstartindices(i):xstartindices(i+1)-1),Ksi,Ksi);
        if nnz(X)<0.05*Ksi2;
            X=sparse(X);
        end
        y(startindices(i):startindices(i+1)-1)=...
            temp'...
            \(X...
            /temp);
    end
else
    y=[];
end