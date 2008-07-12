/*
%                               [Lden,L.d] = dpr1fact(x, d, Lsym, smult, maxu)
% DPR1FACT  Factor d[iag] p[lus] r[ank] 1:
%    [Lden,L.d] = dpr1fact(x, d, Lsym, smult, maxu)
%    Computes fi and d such that
%       diag(d_IN) + x*diag(smult)*x' = 
%(PI_{i=1}^n L(p_OUT^i,beta_i)) * diag(d_OUT) * (PI_{i=1}^n L(p_OUT^i,beta_i))'
%    where L(p,beta) = eye(n) + tril(p*beta',-1).
%    
% Lden.dopiv(k) = 1 if p(:,k) has been reordered, with permutation in
% Lden.pivperm.
% We reorder if otherwise |p(i,k)*beta(j,k)| > maxu.
%
% SEE ALSO fwdpr1,bwdpr1,sedumi
% ******************** INTERNAL FUNCTION OF SEDUMI ********************
function [Lden,L.d] = dpr1fact(x, d, Lsym, smult, maxu)

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

#include <math.h>
#include <string.h>
#include "mex.h"
#include "blksdp.h"

#define LDEN_OUT myplhs[0]
#define D_OUT myplhs[1]
#define NPAROUT 2

#define X_IN prhs[0]
#define D_IN prhs[1]
#define LSYMB_IN prhs[2]
#define SMULT_IN prhs[3]
#define MAXU_IN prhs[4]
#define NPARIN 5

/* ============================================================
   DPR1FACT-subroutines: Compact Cholesky for X = diag(d) + p*p'.
   several versions, to allow sequential or permuted ordering.
   ============================================================ */

/* ************************************************************
   dpr1fact - Compact Cholesky for X = diag(d) + p*p'/t to
       X = L(p,beta) * diag(d_OUT) * L(p,beta)'
       where L(p,beta) = eye(m) + tril(p*beta',-1)
   INPUT:
     n - Order of beta. n = min(m,idep), where idep is the
       1st entry where d(idep) = 0 on input. Caller then needs to finish by
       pivoting on idep by itself.
     mu - mu(m) = 0, mu(i) = max(psqr(i+1:mk)), for i=1:mk-1.
     maxu - Controls stability check: we postpone rows such that
       max(abs(L)) <= maxu.
   UPDATED:
     d    - Length n vector: the diagonal entries. On input, the old ones,
          d(1:n) > 0. On output the updated ones after the factorization.
          Remain positive if t > 0.
     fi   - on input, contains the vector x (=p.^2),
          on output it is such that beta(j) = p(j) / fi(j), for
          j not in ph2psqr.i.
     t - Initial t: set t = 1 for D+p*p', set t = -1 for D-p*p'.
   OUTPUT
     ph2psqr - The postponed rows j, with corresponding psqr(j). Controled
       by maxu.
   REMARK:
     Since L=eye(m)+tril(p*beta'), beta(n-1) and fi(n-1) are useful only
     if m > n: it'll be used in rows n:m-1.
   RETURNS: nph2, number of postponed nodes = length(ph2psqr).
   ************************************************************ */
mwIndex dpr1fact(double *fi, double *d, keydouble *ph2psqr, double *pt, const mwIndex n,
             const double *mu, const double maxu)
{
  mwIndex nph2;
  double dj,fij, muph2, t;
  keydouble p2j;

/* ------------------------------------------------------------
   fi(j) = x(j) + t*d(j),  d_new(j) = fi(j)/t,  tnew = fi(j)/d_old(j)
   Store j in p2j.k
   ------------------------------------------------------------ */ 
  t = *pt;
  nph2 = 0;
  muph2 = 0.0;                 /* muph2 = max(psqr(postponed_nodes)) */
  for(p2j.k = 0; p2j.k < n; p2j.k++){
/* ------------------------------------------------------------
   Step j: remains to factor diag(d(j:end)) + p(j:end)*p(j:end)'/t.
   The pivot is d(j) + p(j)^2/t = (t*d(j)+x(j))/t.
   ------------------------------------------------------------ */
    dj = d[p2j.k];
    p2j.r = fi[p2j.k];                    /* p2j = {j, p_j^2} */
    fij = p2j.r + t*dj;                  /* fi(j) = p_j^2 + t*d_j */
/* ------------------------------------------------------------
   max SQR of below-diag = [pj^2 * max(p(j+1:end).^2)] / t^2
   This should not exceed maxu^2 * pivot^2.
   ------------------------------------------------------------ */
    if(p2j.r * MAX(muph2, mu[p2j.k]) <= SQR(maxu * fij)){
      fi[p2j.k] = fij;                 /* pivot j is stable */
      d[p2j.k]  = fij / t;             /* d(j;NEW) = d_j + (p_j^2 / t). */
      t         = fij / dj;            /* Compute new t for next iter. */
    }
    else{
      ph2psqr[nph2++] = p2j;            /* Postpone to phase 2 */
      muph2 = MAX(muph2, p2j.r);        /* max(ph2psqr.r) */
    }
  }
  *pt = t;
  return nph2;
}


/* ************************************************************
   dpr1factperm - Compact Cholesky for X = diag(d) + p*p' to
       X = L(p,beta) * diag(d_OUT) * L(p,beta)'
       where L(p,beta) = eye(m) + tril(p*beta',-1).
       Follows the sequence given in "perm"; realligns accepted pivots
       from start of "perm", stores rejected ones in ph2psqr.
   INPUT:
     n - Order of beta.  n = min(m,idep), where idep is the
       1st entry where d(idep) = 0 on input. Caller then needs to finish by
       pivoting on idep by itself.
     t - Initial t: set t = 1 for D+p*p', set t = -1 for D-p*p'.
     maxu - Controls stability check: we postpone rows such that
       max(abs(L)) <= maxu.
     mu - max(psqr(perm[i+1:m-1])) for all i=1:n (n <= m). NB: in perm-order.
   UPDATED:
     perm - pivot sequence. Evaluate pivots perm(0:n-1). On output,
       perm(0:n-nph2-1) are the accepted pivots. 
     d    - Length n vector: the diagonal entries. On input, the old ones,
          d(1:n) > 0. On output the updated ones after the factorization.
          Remain positive if t > 0.
     fi   - on input, contains the vector x (=p.^2),
          on output s.t. beta(j) = p(j) / fi(j) for j=perm[0:n-nph2-1].
   OUTPUT
     ph2psqr - The postponed rows j, with corresponding psqr(j). Controled
       by maxu.
   REMARK:
     Since L=eye(m)+tril(p*beta'), beta(n-1) and fi(n-1) are useful only
     if m > n: it'll be used in rows n:m-1.
   RETURNS: nph2, number of postponed nodes = length(ph2psqr).
   ************************************************************ */
mwIndex dpr1factperm(double *fi, double *d, keydouble *ph2psqr, double *pt,
                 mwIndex *perm, const mwIndex n, const double *mu, const double maxu)
{
  mwIndex i, jnz, nph2;
  double dj,fij, muph2, t;
  keydouble p2j;

/* ------------------------------------------------------------
   fi(j) = x(j) + t*d(j),  d_new(j) = fi(j)/t,  tnew = fi(j)/d_old(j)
   Store j in p2j.k
   ------------------------------------------------------------ */ 
  t = *pt;
  nph2 = 0;
  muph2 = 0.0;
  jnz = 0;         /* index into perm_OUT, for accepted pivots */
  for(i = 0; i < n; i++){
    p2j.k = perm[i];
    dj = d[p2j.k];
    p2j.r = fi[p2j.k];                    /* p2j = {j, p_j^2} */
    fij = p2j.r + t*dj;                  /* fi(j) = p_j^2 + t*d_j */
    if(p2j.r * MAX(muph2, mu[i]) <= SQR(maxu * fij)){
      fi[p2j.k] = fij;                   /* pivot j is stable */
      perm[jnz++] = p2j.k;
      d[p2j.k]  = fij / t;             /* d(j;NEW) = d_j + (p_j^2 / t). */
      t         = fij / dj;            /* Compute new t for next iter. */
    }
    else{
      ph2psqr[nph2++] = p2j;            /* Postpone to phase 2 */
      muph2 = MAX(muph2, p2j.r);        /* max(ph2psqr.r) */
    }
  }
  mxAssert(jnz + nph2 == n, "");
  *pt = t;
  return nph2;
}

/* ************************************************************
   ph2dpr1fact - Compact Cholesky for X = diag(d) + p*p' to
       X = L(p,beta) * diag(d_OUT) * L(p,beta)'
       where L(p,beta) = eye(m) + tril(p*beta',-1)
   INPUT:
     n - Order of psqr (number of phase-2 rows).
     t - Initial t: output from 1st phase; is mon. incr.
       t >= 1 for D+p*p', whereas -1 <= t < 0 for D-p*p'.
   UPDATED:
     psqr - Contains the sparse vector (p.^2), where the row-indices
          are the postponed row numbers. On output, the r-values are
          replaced by fi (so that beta = p ./ fi).
     d    - the diagonal entries. On input, the old ones,
          on output the updated ones after the factorization.
          Only those with psqr.i-indices are changed (should be
          all positive already on input).
   REMARK:
     Since L=eye(m)+tril(p*beta'), beta(n-1) and fi(n-1) are useful only
     if m > n: it'll be used in rows n:m-1.
   ************************************************************ */
void ph2dpr1fact(keydouble *psqr, double *d, double *pt, const mwIndex n)
{
  mwIndex j, jnz;
  double dj,fij,t;
  t = *pt;
/* ------------------------------------------------------------
   fi(j) = x(j) + t*d(j),  d_new(j) = fi(j)/t,  tnew = fi(j)/d_old(j)
   ------------------------------------------------------------ */ 
  for(jnz = 0; jnz < n; jnz++){
    j = (psqr+jnz)->k;
    dj = d[j];
    fij = ((psqr+jnz)->r += t*dj);        /* fi(j) = p_j^2 + t*d_j */
    d[j]    = fij / t;             /* d(j;NEW) = d_j + (p_j^2 / t). */
    t       = fij / dj;            /* Compute new t for next iter. */
  }
  *pt = t;
}

/* ============================================================
   MAIN routine for Compact Cholesky for X = diag(d) + p*p'.
   redirects to the dpr1fact subroutines.
   ============================================================ */

/* ************************************************************
   PROCEDURE dodpr1fact - Factors diag +/- rank-1:
     (D+t*p*p')(perm) = L * diag(d_NEW(perm)) * L',
     L = I+tril(p(perm)*beta',-1).
   INPUT
     p    - length m. We've to factor diag(d)+ (1/t) * p*p'.
     t    - scalar: 1 for adding p*p', -1 for subtracting p*p'.
     maxu - scalar >= 1: The factor L(p,beta) = I+tril(p(perm)*beta',-1)
       will be such that max(abs(L)) <= maxu by choosing perm-ordering.
     m    - length(p).
   UPDATED
     d    - length m. The diagonal. This factors
       diag(d_OLD)+t*p*p' = L(p,beta) * diag(d_NEW) * L(p,beta)'
   OUTPUT
     beta - Length <= m (actual length returned in *pm).
     perm - Length m. Only written if RETURN=1, which means that the
       original ordering was not maxu-stable. Pivot ordering on p,d.
     pn   -  *pn = length(beta) <= m; n<m only if there are dependent rows.
     dep  - Length ndep+1. Lists rows i where d(i) == 0. Indices are
       ascending, and dep[ndep] >= m is tail of this list. On output,
       one entry may be removed, and stored in dep[ndep_OLD].
     *pndep - Cardinality of dep. May be decremented on output, if a
       dependency could be removed, i.e. if t > 0 and p(dep) != 0.
   WORK
     psqr - length m float working array, for p.^2 and later "fi".
     kdwork - length m working array for storing postponed
       rows (rowno and psqr(rowno)), which have to be sorted.
   RETURNS 1 if reordered rows into perm; 0 means that we used
     the sequential 0:m-1 ordering.
   CAUTION: If t < 0, one dependency may be added by the
       rank-1 subtraction. The caller should therefore call findnewdep
       afterwards (for t < 0).
   ************************************************************ */
char dodpr1fact(double *beta, mwIndex *perm, double *d, double t, const double *p,
                const mwIndex m, mwIndex *pn, mwIndex *dep, mwIndex *pndep,
                const double maxu, double *psqr, keydouble *kdwork)
{
  mwIndex ndep, n, i, j, nph2, nextj, idep;
  double psqrdep, h;
  double *mu;
  char deldep;

/* ------------------------------------------------------------
   If t = 0, then factor diag(d)+0*p*p' = I*diag(d)*I, i.e. beta=0.
   ------------------------------------------------------------ */
  if(t == 0.0){
    *pn = 0;              /* number of nonzeros in beta */
    return 0;
  }
/* ------------------------------------------------------------
   t is nonzero, replace by tnew := 1/t.
   We've to factor diag(d) + p*p' / tnew.
   ------------------------------------------------------------ */
  t = 1/t;
  ndep = *pndep;
/* ------------------------------------------------------------
   Use beta temporarily as mu(1:m), which lists max(psqr(i+1:m)).
   mu will be used only to select stable pivots, before writing beta.
   ------------------------------------------------------------ */
  mu = beta;
/* ------------------------------------------------------------
   Let psqr = p(1:m).^2
   ------------------------------------------------------------ */
  realHadamard(psqr, p, p, m);
/* ------------------------------------------------------------
   Case A: d(1:mk) > 0 (no dep). Then n = m.
   ------------------------------------------------------------ */
  if(dep[0] >= m){
    *pn = m;
/* ------------------------------------------------------------
   Let mu(m) = 0, mu(i) = max(psqr(i+1:mk)), for i=1:mk-1.
   ------------------------------------------------------------ */
    for(h = 0.0, i = m ; i > 0; i--){
      mu[i-1] = h;
      h = MAX(h, psqr[i-1]);
    }
/* ------------------------------------------------------------
   1st round: pivot sequentially on 1:m, skipping instable ones.
   ------------------------------------------------------------ */
    nph2 = dpr1fact(psqr, d, kdwork, &t, m, mu, maxu);
/* ------------------------------------------------------------
   Write results 1st round: beta = p ./ psqr.
   ------------------------------------------------------------ */
    if(!nph2){                  /* all 1:m handled */
      realHadadiv(beta, p, psqr, m);
      return 0;
    }
    else{                       /* skipped kdwork.k */
      for(i = 0, j = 0; i < nph2; i++){
        nextj = (kdwork+i)->k;
        fromto(perm+j, j, nextj);       /* perm[j-i:nextj-i] = j:nextj */
        realHadadiv(beta + j, p + j, psqr + j, nextj - j);
        j = nextj + 1;              /* skip nextj == (kdwork+i)->k */
        --perm; --beta;             /* keep j valid index */
      }
      fromto(perm+j, j, m);       /* perm[j-i:nextj-i] = j:nextj */
      realHadadiv(beta + j, p + j, psqr + j, m - j);
      perm += m;           /* point just behind accepted pivots */
      beta += m;
/* ------------------------------------------------------------
   Sort rejected nodes in decreasing order of p.^2.
   ------------------------------------------------------------ */
      kdsortdec(kdwork, nph2);
/* ------------------------------------------------------------
   2nd round factorization: ordered.
   ------------------------------------------------------------ */
      ph2dpr1fact(kdwork, d, &t, nph2);
      for(i = 0; i < nph2; i++){
        j = (kdwork+i)->k;
        perm[i] = j;
        beta[i] = p[j] / (kdwork+i)->r;
      }
      return 1;
    } /* if nph2 > 0 */
  } /* if !dep */
/* ------------------------------------------------------------
   If d(1:mk) is NOT positive:
   Let (j,psqrdep) = max{psqr(i) | d(i)==0.0, i=1:m}
   ------------------------------------------------------------ */
  else{
    psqrdep = 0.0;
    for(i = 0; dep[i] < m; i++)
      if(psqr[dep[i]] > psqrdep){
        j = i;
        psqrdep = psqr[dep[i]];
      }
    mxAssert(i <= ndep, "");
/* ------------------------------------------------------------
   Threshold h = maxu^2 * psqrdep
   If all psqr>h have been factorized, we'll pivot on dep[k], if
   t * psqrdep > 0 (otherwise we view this as being zero).
   ------------------------------------------------------------ */
    if(psqrdep > 0.0){        /* we'll remove dependency at idep=dep[j] */
      idep = dep[j];
/* ------------------------------------------------------------
   If psqrdep>0, we can remove dependency idep=dep[j].
   Let dep[j:ndep-1] = dep[j+1:ndep] (incl tail dep[ndep]), then
   let dep[ndep] = idep, and --ndep. For Lorentz cones, removed
   dependencies may get dependent again at the t=-1 step.
   ------------------------------------------------------------ */
      if(t > 0.0){
        deldep = 1;
        memmove(dep+j, dep+j+1, (ndep - j) * sizeof(mwIndex));
        h = SQR(maxu) * psqrdep;
        dep[ndep] = idep;                /* remember removed dependency */
        *pndep = --ndep;
      }
/* ------------------------------------------------------------
   If we're subtracting a rank-1 factor (t<0), then psqrdep should
   be zero (up to rounding errors)
   ------------------------------------------------------------ */
      else{                  /* D - p*p' should be psd, so */
        h = psqrdep;         /* we've to round [0,psqrdep] to 0 */
        deldep = 0;
      }
    }
    else{
      idep = dep[0];           /* psqr(dep) == 0: remains dependent */
      h = 0.0;
      deldep = 0;
    }
/* ------------------------------------------------------------
   PARTITION: perm = [find(psqr > h), idep, remainder].
   Then let n be j = length(find(psqr > h)).
   Temporarily use nph2 = m-length(remainder).
   ------------------------------------------------------------ */
    for(i = 0, j = 0, nph2 = m; i < idep; i++)
      if(psqr[i] > h)
        perm[j++] = i;
      else
        perm[--nph2] = i;
    for(++i; i < m; i++)               /* skip over i = idep */
      if(psqr[i] > h)
        perm[j++] = i;
      else
        perm[--nph2] = i;
    mxAssert(j == nph2-1,"");
    perm[j] = idep;                     /* finally insert idep */
    n = j;                       /* length(find(psqr > h)) */
    *pn = j + deldep;            /* cardinality of beta */
/* ------------------------------------------------------------
   Now h=max(psqr(perm(n+1:m))).
   Let mu(i) = max(psqr(perm(i+1:m))).
   ------------------------------------------------------------ */
    for(i = n ; i > 0; i--){
      mu[i-1] = h;
      h = MAX(h, psqr[perm[i-1]]);
    }
/* ------------------------------------------------------------
   1st round: pivot sequentially on perm(1:n), skipping instable ones.
   The stable pivots are re-alligned at start of perm.
   ------------------------------------------------------------ */
    nph2 = dpr1factperm(psqr, d, kdwork, &t, perm, n, mu, maxu);
/* ------------------------------------------------------------
   Write results 1st round: beta = p(perm(1:n-nph2)) ./ psqr(perm(1:n-nph2)).
   ------------------------------------------------------------ */
    n -= nph2;          /* cardinality 1st round */
    for(i = 0; i < n; i++){
      j = perm[i];
      beta[i] = p[j] / psqr[j];
    }
    perm += n;         /* handled 1st round */
    beta += n;
/* ------------------------------------------------------------
   Sort rejected nodes in decreasing order of p.^2.
   ------------------------------------------------------------ */
    if(nph2){
      kdsortdec(kdwork, nph2);
/* ------------------------------------------------------------
   2nd round factorization: ordered.
   ------------------------------------------------------------ */
      ph2dpr1fact(kdwork, d, &t, nph2);
      for(i = 0; i < nph2; i++){
        j = (kdwork+i)->k;
        perm[i] = j;
        beta[i] = p[j] / (kdwork+i)->r;
      }
    }
/* ------------------------------------------------------------
   If psqrdep > 0, we can now finish off the factorization by
   pivoting on idep == perm[nph2]:
   d_new(i) = p_i^2/t, beta = 1/p_i.
   ------------------------------------------------------------ */
    if(deldep){
      d[idep] = psqr[idep] / t;
      beta[nph2] = 1.0 / p[idep];
    }
  }
  return 1;
}

/* ************************************************************
   PROCEDURE findnewdep - CAUTION: this searches only over previously
     removed dependencies. The rank reduction could however have happened
     elsewehere, viz. last pivot location!!
   INPUT
     ndep    - Number of dependent nodes, d[dep[0:ndep-1]] == 0.
     maxndep - dep is length maxndep+1. dep[ndep+1:maxndep] are previously
          removed dependencies.
     d       - length m vector, m = dep[ndep].
   UPDATED
     dep - length maxndep+1 array. If d[dep[i]] <= 0 for some i > ndep,
       then dep[i] is inserted into dep(0:ndep), so that dep(0:ndep+1) remains
       sorted.
   RETURNS 1 if ndep has to be incremented, i.e. an entry of
     dep(ndep+1:maxndep) is inserted into dep(0:ndep). Otherwise returns 0.
   ************************************************************ */
mwIndex findnewdep(mwIndex *dep, const mwIndex ndep, const mwIndex maxndep, const double *d)
{
  mwIndex i, j, idep;

  for(i = ndep + 1; i <= maxndep; i++)
    if(d[dep[i]] <= 0.0)
      break;
  if(i <= maxndep){
    idep = dep[i];
    j = 0;
    intbsearch(&j, dep, ndep, idep);  /* first j s.t. dep[j] > idep */
    memmove(dep+j+1, dep+j, (i - j) * sizeof(mwIndex));
    dep[j] = idep;
    return 1;
  }
  else
    return 0;
}

/* ============================================================
   PRODFORMFACT does a dpr1fact for each rank-1 update.
   ============================================================ */

/* ************************************************************
   PROCEDURE prodformfact
   INPUT
     xsuper - column k consists of rows 0:xsuper(k+1)-1.
     n      - number of (dense) columns
     smult  - Length n vector. the k-th step adds (D+smult(k)*pk*pk').
     firstpiv - Length n array, first affecting pivot.
     colperm - Length n array, column permutation for smult and firstpiv.
     maxu   - max_k(max abs(Lk)) will be at most maxu. Rows may be
      reordered to achieve this.
   UPDATED
     p  - Length(p) = sum(xsuper). On input, contains the dense columns
       as in X = diag(d) + P*diag(smult(colperm))*P'. On output, a
       product-form forward solve has been made to p(:,2:n).
     d  - length xsuper[n] nonnegative vector. On input, the diagonal w/o dense
       columns. On output, the diagonal in the final product form Cholesky.
     dep - Length ndep+1 list of entries where d(i)=0; dep(0) < dep(1)...;
        dep[ndep] = xsuper[n], the tail.
     pndep  - length of dep, may be decreased on output, if dependencies
       are removed by adding the rank-1 updates..
   OUTPUT
     perm - sum_j(xsuper(j+1)|ordered(j)=1) array, contains a stable pivot
       ordering for those columns where ordered[j]=1.
     beta   - Length length(p). Such that L_k = eye(m) + tril(pk * betak, -1).
     betajc - Length n+1. start of betak. nnz(beta) <= nnz(p).
     ordered - length n. Ordered[j]==1 iff the rows of column j are
       reordered for numerical stability (controled by maxu).
   WORK
     fwork  - length xsuper[n] float working array.
     kdwork - length xsuper[n] (i,r)-working array.
   ************************************************************ */
void prodformfact(double *p, mwIndex *perm, double *beta, mwIndex *betajc,
                  double *d, char *ordered, const mwIndex *xsuper,
                  const mwIndex *colperm, const mwIndex *firstpiv,
                  const double *smult, const mwIndex n, mwIndex *dep, mwIndex *pndep,
                  const double maxu, double *fwork, keydouble *kdwork)
{
  mwIndex k, colk, mk, nk, j, inz, maxndep;
  double *betak, *pk, *pj;
  char useperm;
/* ------------------------------------------------------------
   Initialize. inz points to next avl. place in beta,
   perm is used to store pivot ordering,
   ------------------------------------------------------------ */
  inz = 0;
  maxndep = *pndep;
/* ------------------------------------------------------------
   For all columns k, mk = length(pk), nk = length(betak).
   ------------------------------------------------------------ */
  for(k = 0, pk = p; k < n; k++){
    colk = colperm[k];                     /* pointer into smult, firstpiv */
    betajc[k] = inz;
    mk = xsuper[k+1];
    betak = beta + inz;
    pk += xsuper[k];
    useperm = dodpr1fact(betak, perm, d, smult[colk], pk, mk, &nk, dep, pndep,
                         maxu, fwork, kdwork);
    ordered[k] = useperm;
    if(smult[colk] < 0.0)
      *pndep += findnewdep(dep,*pndep,maxndep,d);
/* ------------------------------------------------------------
   Forward solve on columns p(k+1:n)
   ------------------------------------------------------------ */
    if(smult[colk] != 0.0){
      if(useperm){
        for(j = k+1, pj = pk; j < n; j++){               /* with pivoting */
          pj += xsuper[j];
          if(firstpiv[colperm[j]] <= k)          /*Only if overlapping nzs*/
            fwipr1o(pj, perm, pk, betak, mk, nk);          /* o = ordered */
        }
        perm += mk;    /* full length permutation */
      }
      else
        for(j = k+1, pj = pk; j < n; j++){           /* without pivoting */
          pj += xsuper[j];
          if(firstpiv[colperm[j]] <= k)
            fwipr1(pj, pk, betak, mk, nk);
        }
    }
/* ------------------------------------------------------------
   Point to next column
   ------------------------------------------------------------ */
    inz += nk;
  }
/* ------------------------------------------------------------
   In total, we wrote inz <= length(p) nonzeros in beta.
   ------------------------------------------------------------ */
  betajc[n] = inz;
#ifdef DO_SUPER_SAFE
/* ------------------------------------------------------------
   If smult[i] < 0 for some i, then let dep = find(d<=0), and d(dep) = 0.
   Note: length(d) = m = xsuper[n].
   ------------------------------------------------------------ */
  mk = xsuper[n];
  inz = 0;
  for(j = 0; j < mk; j++)
    if(d[j] <= 0.0){
      d[j] = 0.0;
      dep[inz++] = j;
      mxAssert(inz <= maxndep, "Fatal numerical error in dpr1fact.");
    }
  *pndep = inz;
#endif
}

#define NLDEN_FIELDS 5
/* ============================================================
   MAIN: MEXFUNCTION
   ============================================================ */
/* ************************************************************
   PROCEDURE mexFunction - Entry for Matlab
   ************************************************************ */
void mexFunction(const int nlhs, mxArray *plhs[],
  const int nrhs, const mxArray *prhs[])
{
  mxArray *MY_FIELD;
  mxArray *myplhs[NPAROUT];
  mwIndex m,n,ndep,i,j, permj, pnnz, dznnz, permnnz;
  char *ordered;
  mwIndex *dep, *colperm, *invrowperm, *betajc, *pivperm, *firstpiv;
  double *beta, *d,*betajcPr, *pj, *orderedPr, *fwork, *p, *permPr, *lab;
  const double *colpermPr, *smult, *firstPr;
  const char *LdenFieldnames[] = {"betajc","beta","p","pivperm","dopiv"};
  keydouble *kdwork;
  double maxu;
  jcir x,dz;
/* ------------------------------------------------------------
   Check for proper number of arguments
   ------------------------------------------------------------ */
  mxAssert(nrhs >= NPARIN, "dpr1fact requires more input arguments");
  mxAssert(nlhs <= NPAROUT, "dpr1fact produces less output arguments");
/* ------------------------------------------------------------
   Get inputs (x, lab=d, smult, maxu)
   ------------------------------------------------------------ */
  m = mxGetM(X_IN);                              /* x */
  n = mxGetN(X_IN);
  mxAssert(mxIsSparse(X_IN), "x should be sparse.");
  x.jc = mxGetJc(X_IN);
  x.ir = mxGetIr(X_IN);
  x.pr = mxGetPr(X_IN);
  mxAssert( mxGetM(D_IN) * mxGetN(D_IN) == m, "Size mismatch d.");     /* d */
  mxAssert( mxGetM(SMULT_IN) * mxGetN(SMULT_IN) == n, "Size mismatch smult."); /* smult */
  smult = mxGetPr(SMULT_IN);
  maxu = mxGetScalar(MAXU_IN);                        /* maxu */
/* ------------------------------------------------------------
   DISASSEMBLE structure Lsymb.{dz,perm,first}
   ------------------------------------------------------------ */
  mxAssert(mxIsStruct(LSYMB_IN), "Lsymb should be a structure.");
  MY_FIELD = mxGetField(LSYMB_IN,(mwIndex)0,"dz");      /* Lsymb.dz */
  mxAssert( MY_FIELD != NULL, "Missing field Lsymb.dz.");
  mxAssert(mxGetM(MY_FIELD) == m && mxGetN(MY_FIELD) == n, "Lsymb.dz size mismatch.");
  mxAssert(mxIsSparse(MY_FIELD), "Lsymb.dz must be sparse.");
  dz.jc = mxGetJc(MY_FIELD);
  dz.ir = mxGetIr(MY_FIELD);                             /* (rowperm) */
  MY_FIELD = mxGetField(LSYMB_IN,(mwIndex)0,"perm");        /* Lsymb.perm */
  mxAssert(MY_FIELD != NULL, "Missing field Lsymb.perm.");
  mxAssert(mxGetM(MY_FIELD) * mxGetN(MY_FIELD) == n, "Size mismatch Lsymb.perm."); /* (colperm) */
  colpermPr = mxGetPr(MY_FIELD);
  MY_FIELD = mxGetField(LSYMB_IN,(mwIndex)0,"first");        /* Lsymb.first */
  mxAssert( MY_FIELD != NULL, "Missing field Lsymb.first.");
  mxAssert( mxGetM(MY_FIELD) * mxGetN(MY_FIELD) == n, "Size mismatch Lsymb.first.");
  firstPr = mxGetPr(MY_FIELD);
/* ------------------------------------------------------------
   Let pnnz = sum(dz.jc), dznnz = dz.jc[n].
   ------------------------------------------------------------ */
  for(i = 1, pnnz = 0; i <= n; i++)
    pnnz += dz.jc[i];
  dznnz = dz.jc[n];
/* ------------------------------------------------------------
   Allocate working arrays:
   mwIndex: colperm(n), firstpiv(n), dep(m+1), betajc(n+1), pivperm(pnnz),
     invrowperm(m).
   char: ordered(n)
   double: fwork(dznnz), d(dznnz),
   keydouble: kdwork(dznnz).
   ------------------------------------------------------------ */
  firstpiv= (mwIndex *) mxCalloc(MAX(n,1), sizeof(mwIndex));
  colperm = (mwIndex *) mxCalloc(MAX(n,1), sizeof(mwIndex)); 
  dep     = (mwIndex *) mxCalloc(m+1, sizeof(mwIndex));
  betajc  = (mwIndex *) mxCalloc(n+1, sizeof(mwIndex));
  invrowperm = (mwIndex *) mxCalloc(MAX(m,1),sizeof(mwIndex));
  pivperm = (mwIndex *) mxCalloc(MAX(pnnz,1), sizeof(mwIndex));      /* pivperm */
  ordered = (char *) mxCalloc(MAX(n,1), sizeof(char));       /* boolean */
  fwork   = (double *) mxCalloc(MAX(dznnz,1), sizeof(double));   /* float */
  d = (double *) mxCalloc(MAX(dznnz,1), sizeof(double));
  kdwork  = (keydouble *) mxCalloc(MAX(dznnz,1), sizeof(keydouble)); /*(i,r)*/
/* ------------------------------------------------------------
   ALLOCATE vectors p(pnnz+m), beta(pnnz), .
   NB1: will be assigned to output vectors later.
   NB2: The +m for p is temporary. This will avoid memory problems when
   initializing p(invperm,:) = x, if Lsymb.dz is invalid.
   ------------------------------------------------------------ */
  p = (double *) mxCalloc(MAX(pnnz + m,1), sizeof(double));    /* p */
  beta = (double *) mxCalloc(MAX(pnnz,1), sizeof(double));      /* beta */
/* ------------------------------------------------------------
   Convert colperm and firstpiv  to integer
   ------------------------------------------------------------ */
  for(i = 0; i < n; i++){                 /* colperm(0:n-1) */
    j = colpermPr[i];
    colperm[i] = --j;
  }
  for(i = 0; i < n; i++){
    j = firstPr[i];
    firstpiv[i] = --j;
  }
/* ------------------------------------------------------------
   CREATE  OUTPUT vector lab := dOUT = dIN (duplicate)
   ------------------------------------------------------------ */
  D_OUT = mxDuplicateArray(D_IN);
  lab = mxGetPr(D_OUT);
/* ------------------------------------------------------------
   Let d(1:dznnz) = lab(dz.ir).
   ------------------------------------------------------------ */
  for(i = 0; i < dznnz; i++)
    d[i] = lab[dz.ir[i]];
/* ------------------------------------------------------------
   dep = [find(d<=0), m], ndep = length(find(d==0)
   ------------------------------------------------------------ */
  ndep = 0;
  for(i = 0; i < dznnz; i++)              /* dep = find(d <= 0) */
    if(d[i] <= 0.0)
      dep[ndep++] = i;
  dep[ndep] = m;           /* tail of dep */
/* ------------------------------------------------------------
   Let invrowperm(dz.ir) = 0:dznnz-1, where dznnz = dz.jc[n] <= m
   ------------------------------------------------------------ */
  mxAssert(dznnz <= m,"");
  for(i = 0; i < dznnz; i++)
    invrowperm[dz.ir[i]] = i;
/* ------------------------------------------------------------
   Let p(invrowperm,:) = x(:,colperm)
   ------------------------------------------------------------ */
  for(j = 0, pj = p; j < n; j++){
    pj += dz.jc[j];
    permj = colperm[j];
    for(i = x.jc[permj]; i < x.jc[permj+1]; i++)
      pj[invrowperm[x.ir[i]]] = x.pr[i];
  }
/* ------------------------------------------------------------
   Create output structure Lden
   ------------------------------------------------------------ */
  LDEN_OUT = mxCreateStructMatrix((mwSize)1, (mwSize)1, NLDEN_FIELDS, LdenFieldnames);
/* ------------------------------------------------------------
   Create LDEN.P(pnnz), and realloc p to the size it should have, i.e. pnnz
   ------------------------------------------------------------ */
  MY_FIELD = mxCreateDoubleMatrix(pnnz, (mwSize)1, mxREAL);
  mxSetField(LDEN_OUT, (mwIndex)0,"p", MY_FIELD);
  if(pnnz > 0){
    mxFree(mxGetPr(MY_FIELD));
    if((p = (double *) mxRealloc(p, pnnz * sizeof(double))) == NULL)
      mexErrMsgTxt("Memory allocation error");
    mxSetPr(MY_FIELD, p);
  }
  else
    mxFree(p);
/* ------------------------------------------------------------
   The actual job is done here:
   Adding n rank-1 updates, with a multiple smult(1:n).
   ------------------------------------------------------------ */
  prodformfact(p, pivperm, beta, betajc, d, ordered, dz.jc, colperm,
               firstpiv, smult, n, dep, &ndep, maxu, fwork, kdwork);
/* ------------------------------------------------------------
   THE DIAGONAL IS PERMUTED BACK:
   Bring d back in original ordering: lab(dz.ir) = d(1:dznnz).
   ------------------------------------------------------------ */
  for(i = 0; i < dznnz; i++)
    lab[dz.ir[i]] = d[i];
/* ------------------------------------------------------------
   Let permnnz = sum{dz.jc[j] | ordered[j]==1}, and set
   Lden.pivperm = pivperm (mwIndex to double, but C-form)
   ------------------------------------------------------------ */
  for(i = 0, permnnz = 0; i < n; i++)
    permnnz += ordered[i] * dz.jc[i+1];
  mxAssert(permnnz <= pnnz, "");
  MY_FIELD = mxCreateDoubleMatrix(permnnz, (mwSize)1, mxREAL);
  mxSetField(LDEN_OUT, (mwIndex)0,"pivperm", MY_FIELD);
  permPr = mxGetPr(MY_FIELD);
  for(i = 0; i < permnnz; i++)
    permPr[i] = pivperm[i];                 /* mwIndex to double */
/* ------------------------------------------------------------
   Create LDEN.BETAJC(n+1)
   ------------------------------------------------------------ */
  MY_FIELD = mxCreateDoubleMatrix(n + 1, (mwSize)1, mxREAL);
  mxSetField(LDEN_OUT, (mwIndex)0,"betajc", MY_FIELD);
  betajcPr = mxGetPr(MY_FIELD);
  for(i = 0; i <= n; i++){
    j = betajc[i];
    betajcPr[i] = ++j;
  }
/* ------------------------------------------------------------
   Create LDEN.BETA(betajc[n])
   ------------------------------------------------------------ */
  MY_FIELD = mxCreateDoubleMatrix(betajc[n], (mwSize)1, mxREAL);
  mxSetField(LDEN_OUT, (mwIndex)0,"beta", MY_FIELD);
  if(betajc[n] > 0){
    mxFree(mxGetPr(MY_FIELD));
    if((beta = (double *) mxRealloc(beta, betajc[n] * sizeof(double))) == NULL)
      mexErrMsgTxt("Memory allocation error");
    mxSetPr(MY_FIELD, beta);
  }
  else
    mxFree(beta);
/* ------------------------------------------------------------
   Create LDEN.DOPIV(n)
   ------------------------------------------------------------ */
  MY_FIELD = mxCreateDoubleMatrix(n, (mwSize)1, mxREAL);
  mxSetField(LDEN_OUT, (mwIndex)0,"dopiv", MY_FIELD);
  orderedPr = mxGetPr(MY_FIELD);
  for(i = 0; i < n; i++)
    orderedPr[i] = ordered[i];
/* ------------------------------------------------------------
   Release working arrays
   ------------------------------------------------------------ */
  mxFree(kdwork);
  mxFree(d);
  mxFree(fwork);
  mxFree(ordered);
  mxFree(pivperm);
  mxFree(invrowperm);
  mxFree(betajc);
  mxFree(dep);
  mxFree(colperm);
  mxFree(firstpiv);
/* ------------------------------------------------------------
   Copy requested output parameters (at least 1), release others.
   ------------------------------------------------------------ */
  i = MAX(nlhs, 1);
  memcpy(plhs,myplhs, i * sizeof(mxArray *));
  for(; i < NPAROUT; i++)
    mxDestroyArray(myplhs[i]);
}
