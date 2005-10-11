% FEASREAL Generates a random sparse optimization problem with
%   linear, quadratic and semi-definite constraints. Output
%   can be used by SEDUMI. All data will be real-valued.
%
%   The following two lines are typical:
% > [AT,B,C,K] = FEASREAL;
% > [X,Y,INFO] = SEDUMI(AT,B,C,K);
%
%   An extended version is:
% > [AT,B,C,K] = FEASREAL(m,lpN,lorL,rconeL,sdpL,denfac)
%
% SEE ALSO sedumi AND feascpx.

function [At,b,c,K] = feasreal(m,nLP,qL,rL,nL,denfac)

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
% 02110-1301, USA%

if nargin < 6
    denfac = 0.1;
    if nargin < 5
        nL = 1+ round(7*rand(1,5));
        fprintf('Choosing random SDP-product cone K.s.\n')
        if nargin < 4
            rL = 3 + round(4*rand(1,5));
            fprintf('Choosing random LORENTZ R-product cone K.r.\n')
            if nargin < 3
                qL = 3 + round(4*rand(1,5));
                fprintf('Choosing random LORENTZ-product cone K.q.\n')
                if(nargin < 2)
                    nLP = 10;
                    fprintf('Choosing default K.l=%3.0f.\n',nLP)
                    if nargin < 1
                        m=10;
                        fprintf('Choosing default m=%3.0f.\n',m)
                    end
                end
            end
        end
    end
end
if isempty(nLP)
    nLP = 0;
end
nblk = length(nL) + length(rL) + length(qL) + nLP;
denfac = min(1, max(denfac,2/nblk));
fprintf('Choosing block density denfac=%6.4f per row\n',denfac)
Apattern=sparse(m,nblk);
for j=1:m
    pati = sprandn(1,nblk,denfac);
    Apattern(j,:) = (pati~= 0);
end
sumnLsqr = sum(nL.^2);
x = zeros(nLP+sum(qL)+sum(rL)+sumnLsqr,1);
x(1:nLP) = rand(nLP,1);
At = sparse(length(x),m);
At(1:nLP,:) = sprandn(Apattern(:,1:nLP)');
firstk = nLP+1;
%[nA, mA] = size(At);
% LORENTZ:
for k=1:length(qL)
    nk = qL(k);  lastk = firstk + nk - 1;
    x(firstk) = rand;
    for j=1:m
        if Apattern(j,nLP+k)==1
            At(firstk:lastk,j) = sprand(nk,1,1/sqrt(nk));
        end
    end
    firstk = lastk + 1;
end
% RCONE (Rotated Lorentz):
for k=1:length(rL)
    nk = rL(k);  lastk = firstk + nk - 1;
    x(firstk) = rand; x(firstk+1) = rand;
    for j=1:m
        if Apattern(j,nLP+length(qL)+k)==1
            At(firstk:lastk,j) = sprand(nk,1,1/sqrt(nk));
        end
    end
    firstk = lastk + 1;
end
% SDP:
for k=1:length(nL)
    nk = nL(k);  lastk = firstk + nk*nk - 1;
    Xk = diag(rand(nk,1));            % diagonal, to keep sparsity structure
    x(firstk:lastk) = Xk;             %although symmetric, store nk^2 elts.
    for j=1:m
        if Apattern(j,nLP+length(qL)+length(rL)+k)==1
            Aik = sprandsym(nk,1/nk);       %on average 1 nonzero per row
            At(firstk:lastk,j) = vec( (Aik+Aik') );
        end
    end
    firstk = lastk + 1;
end
b = full(At'*x);
y = rand(m,1)-0.5;
c = sparse(At*y)+x;
K.l = nLP; K.q = qL; K.r = rL; K.s = nL;
