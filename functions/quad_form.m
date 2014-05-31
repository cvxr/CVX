function [ w, success ] = quad_form( x, Q, v, w )

%QUAD_FORM quadratic form.
%   QUAD_FORM(x,Q) is real(x'*Q*x) = x'*((Q+Q')/2)*x.
%   QUAD_FORM(x,Q,v,w) is real(x'*(Q*x+v)+w).
%
%   x must be a row or column vector, and Q must either be a scalar or
%   a square matrix with the same number of rows as x. If supplied, v must
%   be a scalar or a vector of the same size as x, and w must be a scalar.
%   If x is a row vector, then x (and v) are Hermitian transpoed before
%   evaluation.
%  
%   NOTE: The use of QUAD_FORM can often be replaced by a call to NORM. For
%   example, if Q is positive definite, then the constraint
%       quad_form(x,Q) <= 1
%   is equivalent to
%       norm(sqrtm(Q)*x) <= 1
%   Generally speaking, the NORM version will be more reliable and more
%   accurate, so we encourage you to make similar conversions whenever
%   possible. We *strongly* discourage the QP-era practice of converting
%   NORM expressions into quadratic forms. 
%
%   Disciplined convex programming information:
%       QUAD_FORM(x,Q,v,w) is neither convex nor concave in x and (Q,v)
%       jointly, so at least one of the two must be constant.
%
%       If (Q,v) is constant, then QUAD_FORM is convex if Q is positive
%       semidefinite, and concave if Q is negative semidefinite. An error 
%       is generated if Q is indefinite (unless x is also constant). 
%       QUAD_FORM is nonmonotonic in x, so x must be affine.
%       
%       If x is constant, then QUAD_FORM is affine in Q, v, and w. The
%       signs of x will govern whether the elements of Q, v, and w may
%       be convex, concave, or affine.

tol = 16 * eps;
tolLDL = 4 * eps;
if nargin < 4,
    w = 0;
    if nargin < 3,
        v = 0;
    end
end

sx = size( x );
if length( sx ) ~= 2 || all( sx ~= 1 ),
    cvx_throw( 'The first argument must be a vector.' );
end
nx = prod( sx );

sQ = size( Q );
if length( sQ ) ~= 2 || sQ( 1 ) ~= sQ( 2 ),
    cvx_throw( 'The second argument must be a scalar or a square matrix.' );
elseif sQ( 1 ) ~= 1 && sQ( 1 ) ~= nx,
    cvx_throw( 'The size of Q is incompatible with the size of x.' );
end

if nargin < 3,
    v = sparse( sx(1), sx(2) );
elseif ~isequal( size( v ), sx ),
    cvx_throw( 'The size of v is incompatible with the size of x.' );
end

if nargin < 4,
    w = 0;
elseif numel(w) ~= 1 || ~isreal(w),
    cvx_throw( 'The fourth argument must be a real scalar.' );
end

if sx(1) ~= nx,
    x = x';
    v = v';
end

success = true;
cvx_optval = [];
if cvx_isconstant( x ),
    x = cvx_constant( x );
    cvx_optval = real( x' * Q * x ) + sum( real( v' * x ) ) + w;
    return
elseif ~cvx_isconstant( Q ) || ~cvx_isconstant( v ),
    cvx_throw( 'Either x or (Q,v) must be constant.' );
elseif ~cvx_isaffine( x ),
    cvx_throw( 'First argument must be affine.' );
end

Q = cvx_constant( Q );
Q = 0.5 * ( Q + Q' );
v = cvx_constant( v );
dQ = diag( Q );

if nnz( Q ) == 0,
    cvx_optval = real( v' * x ) + w;
    return
end

if sQ( 1 ) == 1,
    x = x + 0.5 * ( v / Q );
    w = w - 0.25 * ( v' * v ) / Q;
    if ~isreal( x ), x = abs( x ); end
    w = real( Q ) * sum_square( x ) + w;
    return
end

while true,
    
    %
    % Remove zero rows and columns from Q. If a diagonal element of Q is
    % zero but there are elements on that row or column that are not,
    % then we know that neither Q nor -Q is PSD.
    %
    
    dQ = diag( Q );
    trQ = sum( dQ );
    if ~all( dQ ),
        nnzQ = nnz( Q );
        tt = dQ ~= 0;
        Q = Q( tt, tt );
        if nnz( Q ) ~= nnzQ,
            success = false;
            break
        end
        dQ = dQ( tt );
        w = w + real( v( ~tt, : )' * cvx_subsref( x, ~tt, ':' ) );
        v = v( tt, : );
        x = cvx_subsref( x, tt, ':' );
        sx = length( x );
    end
    
    %
    % Determine the sign of the elements of Q. If they are not all of
    % the same sign, then neither Q nor -Q is PSD. Note that trQ has
    % preserved the sign of our quadratic form, so setting Q=-Q here
    % in the concave case does not cause a problem.
    %

    dQ = dQ > 0;
    if ~all( dQ ),
        if any( dQ ),
            success = false;
        else
            Q = -Q;
        end
    end
    
    %
    % We've had to modify this portion of the code because MATLAB has
    % removed support for the CHOLINC function.
    %
    % First, try a Cholesky. If it successfully completes its
    % factorization without fail, we accept it without question. If
    % it terminates early, we perform a numerical test to see if the
    % result still approximates the square root to good precision.
    %
    % If Cholesky fails, then we assume the matrix is either rank
    % deficient or indefinite. For sparse matrices, we perform an LDL
    % factorization, and remove the contributions of any 2x2 blocks,
    % negative 1x1 blocks, and near-zero 1x1 blocks on the diagonal.
    % If there are no such blocks, we accept it without question; if
    % so, we perform the same numerical test. If the test fails, we 
    % assume, for sparse matrices, at least, that the matrix is
    % indefinite.
    %
    % If the matrix is dense, our final test is an eigenvalue
    % decomposition, the most expensive but the most accurate.
    %

    spQ = nnz(Q) <= 0.1 * nx * nx;
    if spQ,
        Q = sparse( Q );
        [ R, p, prm ] = chol( Q, 'upper', 'vector' );
        if any( diff(prm) ~= 1 ),
            R( :, prm ) = R; %#ok
        end
    else
        Q = full( Q );
        [ R, p ] = chol( Q, 'upper' );
        if p > 1, 
            R = [ R , R' \ Q(1:p-1,p:end) ]; %#ok
        end
    end
    valid = p == 0;
    if ~valid,
        tolQ = tol * norm( Q, 'fro' );
        if p > 1,
            valid = norm( Q - R' * R, 'fro' ) < tolQ;
        end
    end
    if ~valid && spQ,
        [ R, DD, prm ] = ldl( sparse( Q ), 'upper', 'vector' );
        if nnz( R ) > max( sx, 0.1 * sx * ( sx + 1 ) / 2 ), spQ = false; end
        tt = diag(DD,1) == 0;
        tt = [ tt ; true ] & [ true ; tt ] & ( diag(DD) > tolLDL * trQ );
        DD = diag(DD);
        R  = bsxfun( @times, sqrt( DD(tt,:) ), R(tt,:) );
        if any( diff(prm) ~= 1 ), R( :, prm ) = R; end
        valid = all( tt ) || norm( Q - R' * R, 'fro' ) < tolQ;
    end
    if ~valid && ~spQ,
        [ V, D ] = eig( full( Q ) );
        if nnz( V ) <= max( length(V), 0.1 * numel(V) ), V = sparse(V); end
        D = diag( D );
        if any( D(2:end) < D(1:end-1) ),
            [D,ndxs] = sort(D);
            V = V(:,ndxs);
        end
        valid = D(1) > -tol * D(end);
        if valid,
            nzero = nnz( cumsum(D) < tol * abs(trQ) );
            V = V(:,nzero+1:end);
            D = sqrt(D(nzero+1:end));
            R = diag(sparse(D)) * V';
        end
    end
    if ~valid,
        success = false;
        break;
    end
    
    %
    % Scale so that the mean eigenvalue of (1/alpha)*R'*R is one. 
    % Hopefully this will minimize scaling issues.
    %
   
    alpha = trQ / size(R,1);
    w = w + alpha * sum_square_abs( ( R * x ) / sqrt(alpha) ) + real( v' * x );
    break;
    
end

if ~success && nargout == 1,
    cvx_throw( 'Disciplined convex programming error:\n    The second argument must be positive or negative semidefinite.' );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
