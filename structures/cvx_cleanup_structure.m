function xi = cvx_cleanup_structure( xi )

%CVX_CLEANUP_STRUCTURE   Canonicalizes a matrix structure basis.
%    CVX_CLEANUP_STRUCTURE(X), where X is an m x n matrix, converts X
%    to row reduced echelon form with the zero rows removed. The matrices
%    X are assumed to come from CVX's matrix structure facility, and as
%    such certain assumptions are made about both X and its RREF: that it
%    is sparse and its nonzero elements are ratios of small integers. As a
%    result this routine should not be used for general RREF computations
%    without modifying it to remove some of the cleanup code.

% Reduce using an LU factorization
[LL,xi,PP] = lu(xi); %#ok
[m,n] = size(xi);

% Remove the entries that are close to zero
tol = 16 * eps;
xi  = xi .* ( abs(xi) > tol * norm(xi,'inf') );

% Find the locations of the leading element in each row. To do this we first
% find the first element in each row. Transposing xi insures that the
% indices are sorted properly to accomplish this.
[jj,ii] = find(xi');
if isempty(jj),
    xi = sparse(0,n);
    return
end
dd = [true;diff(ii)~=0];
ii = ii(dd);
jj = jj(dd);

% Sort the rows so that the leftmost nonzero is first (the LU factorization
% does this already much of the time, but in rank-degenerate cases further
% sorting is needed.) From that select a unique set of columns to use.
[jj,jndx] = sort(jj);
dd = [true;diff(jj)~=0];
ii = ii(jndx(dd));
jj = jj(dd);

% Divide through the rows ii by this full-rank triangle. The other rows
% must already be zero or be dependent upon these rows, so we remove them.
% Q is the result except that its columns are scrambled.
rr = length(ii);
j2 = (1:n)'; j2(jj) = [];
Q  = xi(ii,jj)\xi(ii,j2);

% Reduce roundoff error by converting the values to ratios of integers.
[i3,j3,vv] = find(Q);
[vn,vd] = rat(vv,tol);
xi = sparse([(1:rr)';i3],[jj;j2(j3)],[ones(rr,1);vn./vd],rr,n);

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
