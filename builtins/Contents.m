% CVX: Built-in operators and functions supported in CVX models.
%   
%    The following operators and functions are included with MATLAB but
%    have been extended to be used in CVX models as well. Typing
%       help cvx/<func>
%    where <func> is one of the names listed below will provide specific
%    help on the proper use of that item in CVX models---including any
%    restrictions imposed by the DCP and DGP rulesets.
%
% Computational operators:
%    plus/uplus (+), minus/uminus (-), times (.*), mtimes (*), 
%    ldivide (.\), mldivide (\), rdivide (./), mrdivide (/), 
%    power (.^), mpower (^), subsref/subsasgn/end (())
% Relational operators:
%    eq (==), ge (>=), gt (>), le (<=), lt(<), ne (~=).
% Linear/affine functions:
%    blkdiag, cat, conj, ctranspose, cumsum, diag, find, hankel, horzcat,
%    imag, kron, permute, polyval, real, reshape, sparse, sum, toeplitz,
%    tril, triu, vertcat
% Nonlinear functions:
%    abs, max, min, norm, prod, sqrt
% Query functions:
%    disp/display, end, isempty, isequal, length, isreal, ndims, nnz,
%    numel, size, spy
