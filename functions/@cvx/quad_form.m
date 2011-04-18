function [ cvx_optval, success ] = quad_form( x, Q, v, w )

%QUAD_FORM   Internal cvx version.

%
% Check sizes and types
%

error( nargchk( 2, 4, nargin ) );
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

    if sv < sx, v = v(ones(sx,1),1); end
    cvx_optval = 0;
    while true,
        
        %
        % Remove zero rows and columns from Q. If a diagonal element of Q is
        % zero but there are elements on that row or column that are not,
        % then we know that neither Q nor -Q is PSD.
        %
        
        Q = 0.5 * ( Q + Q' );
        dQ = diag( Q );
        trQ = sum(dQ);
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
                cvx_optval = real( v( ~tt, : )' * cvx_subsref( x, ~tt, ':' ) );
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
        % First, try a Cholesky. If rank deficiency is detected, we may
        % be able to recover a valid square root nonetheless. So we'll try,
        % and if it is accurate to a tight tolerance, we'll use it.
        %

        vbar = 0;
        valid = false;
        if cvx_use_sparse( Q ),
            Q = sparse( Q );
            prm = symamd( Q );
            R = cholinc( Q( prm, prm ), 'inf' );
            R( :, prm ) = R;
            tt = any( isinf( R ), 2 );
            valid = ~any( tt );
            if ~valid, 
                R( tt, : ) = []; 
            end
        else
            Q = full( Q );
            [ R, p ] = chol( Q );
            valid = p == 0;
            if ~valid,
                R = [ R , R' \ Q(1:p-1,p:end) ];
            end
        end
        if ~valid,
            valid = normest( Q - R' * R ) < tol * normest( Q );
        end
        if valid && nnz(v),
            vbar = R' \ v;
        end
        
        %
        % If the Cholesky fails, use an eigenvalue decompositon.
        %
        
        if ~valid,
            [ V, D ] = eig( full( Q ) );
            if cvx_use_sparse( V ), 
                V = sparse( V ); 
            end
            D = diag( D );
            Derr = tol * max( D );
            if min( D ) < - Derr, 
                success = false;
                break; 
            end
            tt = D > Derr;
            V = V( :, tt );
            D = D( tt );
            R = diag(sparse(D)) * V';
            if nnz(v),
                vbar = D .\ ( V' * v );
            end
        end
        
        %
        % Scale so that the mean eigenvalue of (1/alpha)*R'*R is one. 
        % Hopefully this will minimize scaling issues.
        %
       
        alpha = trQ / size(R,1);
        if nnz(v),
            v = v - R' * vbar;
            vbar = 0.5 * sg * vbar;
        end
        wbar = w - sg * vbar' * vbar;
        cvx_optval = cvx_optval + sg * alpha * sum_square_abs( ( R * x + vbar ) / sqrt(alpha) ) + real( v' * x ) + wbar;
        break;
        
    end
    
    if ~success,
        error( 'The second argument must be positive or negative semidefinite.' );
    end

end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
