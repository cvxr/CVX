/*
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

*/
#include "blksdp.h"

/* ------------------------------------------------------------
   PROCEDURE fwipr1 - I(dentity) P(lus) R(ank)1 forward solve.
   INPUT:
   p    - length m floats
   beta - length n floats
   m, n - order of p and beta, resp.
   UPDATED:
   y - Length m. On input, contains the rhs. On output, the solution to
         L(p,beta)*yNEW = yOLD
   ------------------------------------------------------------ */
void fwipr1(double *y, const double *p, const double *beta,
            const mwIndex m, const mwIndex n)
{
  mwIndex i;
  double yi,betai,t;

  if(n < 1)           /* If L = I, y remains the same */
    return;
/* ------------------------------------------------------------
   Solve (eye(m) + tril(p*beta',-1)) * yNEW = yOLD,
   where beta[n:m-1] = 0.
   ------------------------------------------------------------ */
  yi = y[0];
  betai = beta[0];
  for(t=0.0, i=1; i < n; i++){
/* ------------------------------------------------------------
   Let t = beta(0:i-1)'*y(0:i-1). Then solve yi from
   t * p(i) + y(i) = yOLD(i)
   ------------------------------------------------------------ */
    t += yi * betai;
    yi = (y[i] -= t * p[i]);
    betai = beta[i];
  }
  if(n < m){
    t += yi * betai;
/* ------------------------------------------------------------
   For i=n:m-1, t remains unchanged.
   ------------------------------------------------------------ */
    addscalarmul(y+i, -t, p+i, m-n);
  }
}

/* ------------------------------------------------------------
   PROCEDURE fwipr1o - I(dentity) P(lus) R(ank)1 forward solve, O(rdered).
   INPUT:
   perm   - length m permutation on p and y.
   p       - length m floats, corresponding to original indices (unpermuted).
   beta - length n floats, corresponding to indices in pir (i.e. already
         permuted); n <= m.
   m, n - order of p and beta, resp.; n <= m.
   UPDATED:
   y - Length m. On input, contains the rhs. On output, the solution to
         L(p(perm),beta)*yNEW(perm) = yOLD(perm)
   ------------------------------------------------------------ */
void fwipr1o(double *y, const mwIndex *perm, const double *p, const double *beta,
             const mwIndex m, const mwIndex n)
{
  mwIndex i, permi;
  double yi,betai,t;

  if(n < 1)           /* If L = I, y remains the same */
    return;
/* ------------------------------------------------------------
   Solve (eye(m) + tril(p*beta',-1)) * yNEW(perm) = yOLD(perm),
   where beta[n:m-1] = 0.
   ------------------------------------------------------------ */
  yi = y[perm[0]];
  betai = beta[0];
  for(t=0.0, i=1; i < n; i++){
/* ------------------------------------------------------------
   Let t = beta(0:i-1)'*y(perm(0:i-1)). Then solve yi from
   t * p(perm(i)) + y(perm(i)) = yOLD(perm(i))
   ------------------------------------------------------------ */
    t += yi * betai;
    permi = perm[i];
    yi = (y[permi] -= t * p[permi]);
    betai = beta[i];
  }
  if(n < m){
    t += yi * betai;
/* ------------------------------------------------------------
   For i=n:m-1, t remains unchanged.
   ------------------------------------------------------------ */
    for(; i < m; i++){
      permi = perm[i];
      y[permi] -= t * p[permi];
    }
  }
}
