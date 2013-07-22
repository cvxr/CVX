%                                           y = psdscale(ud,x,K [,transp])
% PSDSCALE  Computes length lenud (=sum(K.s.^2)) vector y.
%   !transp (default) then y[k] = vec(Ldk' * Xk * Ldk)
%   transp == 1 then y[k] = vec(Udk' * Xk * Udk)
%   Uses pivot ordering ud.perm if available and nonempty.
%
% **********  INTERNAL FUNCTION OF SEDUMI **********
%
% See also scaleK, factorK.

function y = psdscale(ud,x,K,varargin)

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

%The function is quite dirty as the type and dimension of the inputs may
%change, that is why we have so many subcases.

if ~isempty(varargin)
    transp=varargin{1};
else
    transp=0;
end

Ks=K.s;
if ~isempty(Ks)
    N=sum(Ks.^2);
    y(N,1)=0;
else
    y=[];
    return
end
startindices=K.sblkstart-K.mainblks(end)+1;
%Sometimes x containts only the PSD part, sometimes the whole thing
xstartindices=startindices+(length(x)-N);

permindices=cumsum([1,Ks]);
if ~transp
    if isstruct(ud)
        if isempty(ud.perm)
            for i=1:K.rsdpN
                Ksi=Ks(i);
                temp=tril(reshape(ud.u(startindices(i):startindices(i+1)-1),Ksi,Ksi));
                y(startindices(i):startindices(i+1)-1)=...
                    temp'...
                    *reshape(x(xstartindices(i):xstartindices(i+1)-1),Ksi,Ksi)...
                    *temp;
            end
        else
            for i=1:K.rsdpN
                Ksi=Ks(i);
                temp=tril(reshape(ud.u(startindices(i):startindices(i+1)-1),Ksi,Ksi));
                X=reshape(x(xstartindices(i):xstartindices(i+1)-1),Ksi,Ksi);
                %We only permute if we have to
                if isequal(ud.perm(permindices(i):permindices(i+1)-1),(1:Ksi)')
                    if nnz(X)/(Ksi*Ksi)<0.1
                        y(startindices(i):startindices(i+1)-1)=...
                            temp'...
                            *sparse(X)...
                            *temp;
                    else
                        y(startindices(i):startindices(i+1)-1)=...
                            temp'...
                            *X...
                            *temp;
                    end
                else
                    y(startindices(i):startindices(i+1)-1)=...
                        temp'...
                        *X(ud.perm(permindices(i):permindices(i+1)-1),ud.perm(permindices(i):permindices(i+1)-1))...
                        *temp;
                end
            end
        end
    else
        for i=1:K.rsdpN
            Ksi=Ks(i);
            temp=tril(reshape(ud(startindices(i):startindices(i+1)-1),Ksi,Ksi));
            X=reshape(x(xstartindices(i):xstartindices(i+1)-1),Ksi,Ksi);
            if spars(X)<0.1
                y(startindices(i):startindices(i+1)-1)=...
                    temp'...
                    *sparse(X)...
                    *temp;
            else
                y(startindices(i):startindices(i+1)-1)=...
                    temp'...
                    *X...
                    *temp;
            end
        end
    end
else %transp
    if isstruct(ud)
        if isempty(ud.perm)
            for i=1:K.rsdpN
                Ksi=Ks(i);
                temp=triu(reshape(ud.u(startindices(i):startindices(i+1)-1),Ksi,Ksi));
                y(startindices(i):startindices(i+1)-1)=...
                    temp'...
                    *reshape(x(xstartindices(i):xstartindices(i+1)-1),Ksi,Ksi)...
                    *temp;
            end
        else
            for i=1:K.rsdpN
                Ksi=Ks(i);
                temp=triu(reshape(ud.u(startindices(i):startindices(i+1)-1),Ksi,Ksi));
                X=reshape(x(xstartindices(i):xstartindices(i+1)-1),Ksi,Ksi);
                if isequal(ud.perm(permindices(i):permindices(i+1)-1),(1:Ksi)')
                    if nnz(X)/(Ksi*Ksi)<0.1
                        y(startindices(i):startindices(i+1)-1)=temp'...
                            *sparse(X)...
                            *temp;
                    else
                        y(startindices(i):startindices(i+1)-1)=temp'...
                            *X...
                            *temp;
                    end
                else
                    tempy=zeros(Ksi,Ksi);
                    tempy(ud.perm(permindices(i):permindices(i+1)-1),ud.perm(permindices(i):permindices(i+1)-1))=temp'...
                        *X...
                        *temp;
                    y(startindices(i):startindices(i+1)-1)=tempy;
                end
            end
        end
    else
        for i=1:K.rsdpN
            Ksi=Ks(i);
            temp=triu(reshape(ud(startindices(i):startindices(i+1)-1),Ksi,Ksi));
            y(startindices(i):startindices(i+1)-1)=...
                temp'...
                *reshape(x(xstartindices(i):xstartindices(i+1)-1),Ksi,Ksi)...
                *temp;
        end
    end
end