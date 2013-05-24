function [ cvx_optval, success ] = quad_form( x, Q, v, w )

%QUAD_FORM   Internal cvx version.

%
% Check sizes and types
%

error( nargchk( 2, 4, nargin ) ); %#ok
tol = 10 * eps;
if nargin < 4,
    w = 0;
    if nargin < 3,
        v = 0;
    end
end
sx = size( x );
if length( sx ) ~= 2 || all( sx ~= 1 ),
    error( 'The first argument must be a vector.' );
else
    sx = prod( sx );
end
sQ = size( Q );
if length( sQ ) ~= 2 || sQ( 1 ) ~= sQ( 2 ),
    error( 'The second argument must be a scalar or a square matrix.' );
elseif all( sQ( 1 ) ~= [ 1, sx ] ),
    error( 'The size of Q is incompatible with the size of x.' );
end
sv = size( v );
if length( sv ) > 2 || all( sv ~= 1 ),
    error( 'The third argument must be a vector.' );
elseif all( prod( sv ) ~= [ 1, sx ] ),   
    error( 'The size of v is incompatible with the size of x.' );
else
    sv = prod( sv );
end
if numel( w ) ~= 1,
    error( 'The fourth argument must be a real scalar.' );
end

w = real( w );
v = vec( v );
x = vec( x );
success = true;
if cvx_isconstant( x ),
    
    %
    % Constant x, affine Q
    %

    x = cvx_constant( x );
    cvx_optval = real( x' * Q * x ) + sum( real( v' * x ) ) + w;
    return

elseif ~cvx_isaffine( x ),

    error( 'First argument must be affine.' );
    
elseif ~cvx_isconstant( Q ) || ~cvx_isconstant( v ),
    
    error( 'Either x or (Q,v) must be constant.' );
    
end
Q = cvx_constant( Q );
v = cvx_constant( v );
if nnz( Q ) == 0,
    
    %
    % Zero Q, affine x
    %

    cvx_optval = sum( real( v' * x ) ) + w;
    return
    
elseif sQ( 1 ) == 1,
    
    %
    % Constant scalar Q, affine x
    %
    
    x = x + 0.5 * ( v / Q );
    w = w - 0.25 * ( v' * v ) / Q;
    cvx_optval = real( Q ) * sum_square_abs( x ) + w;
    return
    
else

    %
    % Constant matrix Q, affine x
    %

    if sv < sx, 
        v = v(ones(sx,1),1); 
    end
    cvx_optval = w;
    while true,
        
        %
        % Remove zero rows and columns from Q. If a diagonal element of Q is
        % zero but there are elements on that row or column that are not,
        % then we know that neither Q nor -Q is PSD.
        %
        
        Q = 0.5 * ( Q + Q' );
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
            if nnz( v ),
                cvx_optval = cvx_optval + real( v( ~tt, : )' * cvx_subsref( x, ~tt, ':' ) );
            end
            v = v( tt, : );
            x = cvx_subsref( x, tt, ':' );
        end
        
        %
        % Determine the sign of the elements of Q. If they are not all of
        % the same sign, then neither Q nor -Q is PSD.
        %

        dQ = dQ > 0;
        if all( dQ ),
            sg = +1;
        elseif any( dQ ),
            success = false;
            break
        else
            sg = -1;
            Q = -Q;
        end
        
        %
        % We've had to modify this portion of the code because MATLAB has
        % removed support for the CHOLINC function.
        %
        % First, try a Cholesky. If it successfully completes its
        % factorization without fail, we will assume that it is valid
        % with no further checking.
        %
        % If Cholesky fails, then we have detected either rank deficiency
        % or an indefinite matrix. In the non-sparse case, we attempt to
        % construct a square root from the halted Cholesky results. In the
        % sparse case, we perform an LDL, deleting any 2x2 blocks and any
        % nonpositive 1x1 blocks.
        %
        % Neither of these approaches is guaranteed to produce a valid
        % square root in the semidefinite case. So we must test the result
        % and reject any failures.
        %
        % We have had to remove the incomplete Cholesky factorization
        % approach because MATLAB has stopped supporting it. We have
        % replaced it with an LDL. We delete any 2x2 blocks and any
        % nonpositive 1x1 blocks. We do not guarantee this will always
        % produce a numerically valid square root in the semidefinite case,
        % but our norm test can reject failures.
        %

        if cvx_use_sparse( Q ),
            Q = sparse( Q );
            [ R, p, prm ] = chol( Q, 'upper', 'vector' );
            if p ~= 0,
                [ R, DD, prm ] = ldl( Q, 'upper', 'vector' );
                tt = diag(DD,1) == 0;
                tt = [ tt ; true ] & [ true ; tt ] & diag(DD) > 0;
                DD = diag(DD);
                R  = bsxfun( @times, sqrt(DD(tt,:)), R(tt,:) );
            end
            R( :, prm ) = R;
        else
            [ R, p ] = chol( full( Q ), 'upper' );
            if p ~= 0,
                R = [ R , R' \ Q(1:p-1,p:end) ];
            end
        end
        valid = p == 0;
        if ~valid,
            valid = norm( Q - R' * R, 'fro' ) < tol * norm( Q, 'fro' );
        end
        
        %
        % If Cholesky and LDL fail, do an eigendecomposition.
        %
        
        if ~valid,
            [ V, D ] = eig( full( Q ) );
            if cvx_use_sparse( V ), 
                V = sparse( V ); 
            end
            D = diag( D );
            if any( D(2:end) < D(1:end-1) ),
                [D,ndxs] = sort(D);
                V = V(:,ndxs);
            end
            if D(1) < -tol * D(end),
                success = false;
                break;
            end
            nzero = nnz(cumsum(D)<tol*abs(trQ));
            V = V(:,nzero+1:end);
            D = sqrt(D(nzero+1:end));
            R = diag(sparse(D)) * V';
        end
        
        %
        % Scale so that the mean eigenvalue of (1/alpha)*R'*R is one. 
        % Hopefully this will minimize scaling issues.
        %
       
        alpha = trQ / size(R,1);
        cvx_optval = cvx_optval + alpha * sum_square_abs( ( R * x ) / sqrt(alpha) ) + real( v' * x );
        break;
        
    end
    
    if ~success,
        error( 'The second argument must be positive or negative semidefinite.' );
    end

end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
