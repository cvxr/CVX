%                                                     [lab,q,f] = eigK(x,K)
%
% EIGK  Computes the spectral values ("eigenvalues") or even the complete
%       spectral decomposition  of a vector x with respect to a self-dual
%       homogeneous cone K.
%
% > LAB = EIGK(x,K) This yield the spectral values of x with respect to
%       the self-dual homogeneous cone that you describe in the structure
%       K. Up to 3 fields can be used, called K.l, K.q and K.s, for
%       Linear, Quadratic and Semi-definite. Type `help sedumi' for more
%       details on this structure.
%
%       The length of the vector LAB is the order of the cone. Remark that
%       x in K if and only if LAB>=0, and x in int(K) if and only if LAB>0.
%
% > [LAB,Q,F] = EIGK(x,K) Also produces the eigenvectors for the symmetric and
%       Hermitian parts of x (corresponding to K.s), and the Lorentz frame F
%       for the Lorentz blocks in x (corresponding to K.q). This extended
%       version of EIGK is intended for INTERNAL USE BY SEDUMI.
%
% See also sedumi, mat, vec, eyeK.

function [lab,q,f] = eigK(x,K) %#ok

% THE M-FILE VERSION OF THIS FUNCTION IS HERE ONLY AS ILLUSTRATION.
% SEE THE C-SOURCE FOR THE MEX-VERSION.

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

disp('The SeDuMi binaries are not installed.')
disp('In Matlab, launch "install_sedumi" in the folder you put the SeDuMi files.')
disp('For more information see the file Install.txt.')
error(' ')

lab = zeros(K.l + 2*length(K.q) + sum(K.s),1); %#ok
% ----------------------------------------
% LP: lab = x
% ----------------------------------------
lab(1:K.l) = x(1:K.l);
firstk = K.l+1;
nextlab = K.l+1;
% ----------------------------------------
% LORENTZ: (lab, f) = qeig(x)
% ----------------------------------------
if nargout < 3
    for k=1:length(K.q)
        nk = K.q(k); lastk = firstk + nk - 1;
        lab(nextlab:nextlab+1) = qeig(x(firstk:lastk));
        firstk = lastk + 1; nextlab = nextlab + 2;
    end
else
    nextf = 1;
    for k=1:length(K.q)
        nk = K.q(k); lastk = firstk + nk - 1;
        [labk, fk] = qeig(x(firstk:lastk));
        lab(nextlab:nextlab+1) = labk;
        f(nextf:nextf+nk-2) = fk;
        firstk = lastk + 1; nextlab = nextlab + 2; nextf = nextf + nk-1;
    end
end
% ----------------------------------------
% SDP: [q,lab] = eig(x)
% ----------------------------------------
offq = firstk - 1;
q = zeros(length(x)-offq,1);
for k = 1:K.rsdpN
    nk = K.s(k); lastk = firstk + nk*nk - 1;
    Xk = mat(x(firstk:lastk),nk);
    Xk = (Xk + Xk') / 2;
    [Qk, Labk] = eig(Xk);
    lab(nextlab:nextlab+nk-1) = diag(Labk);
    q(firstk-offq:lastk-offq) = Qk;
    firstk = lastk + 1; nextlab = nextlab + nk;
end
for k = (1+K.rsdpN):length(K.s)
    nk = K.s(k); ifirstk = firstk + nk*nk; lastk = ifirstk + nk*nk - 1;
    Xk = mat(x(firstk:ifirstk-1)+sqrt(-1)*x(ifirstk:lastk),nk);
    [Qk, Labk] = eig(Xk);
    lab(nextlab:nextlab+nk-1) = real(diag(Labk));
    q(firstk-offq:ifirstk-1-offq) = real(Qk);
    q(ifirstk-offq:lastk-offq) = imag(Qk);
    firstk = lastk + 1; nextlab = nextlab + nk;
end