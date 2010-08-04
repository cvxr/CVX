function [ cvx_optval, success ] = quad_form( x, Q, v, w )

%QUAD_FORM   Internal cvx version.

%
% Check sizes and types
%

error( nargchk( 2, 4, nargin ) );
tol = 4 * eps;
if nargin < 4,
    w = 0;
    if nargin < 3,
        v = 0;
    end
end
sx = size( x );
if length( sx ) ~= 2 || all( sx > 1 ),
    error( 'The first argument must be a row or column.' );
else
    sx = prod( sx );
end
sQ = size( Q );
if length( sQ ) ~= 2 || sQ( 1 ) ~= sQ( 2 ),
    error( 'The second argument must be a scalar or a square matrix.' );
elseif sQ( 1 ) ~= sx && sQ( 1 ) ~= 1,
    error( 'Sizes are incompatible.' );
end
sv = size( v );
if all( prod( sv ) ~= [ 1, sx ]  ),
    error( 'Sizes are incompatible.' );
elseif nnz( sv ~= 1 ) > 1,
    error( 'The third argument must be a vector.' );
end
if numel( w ) ~= 1,
    error( 'The fourth argument must be a scalar.' );
end

v = vec( v );
x = vec( x );
success = true;
if cvx_isconstant( x ),

    if isreal( Q ) || isreal( x ),

        %
        % Constant x, affine Q, real case
        %

        x = real( x );
        Q = real( Q );
        cvx_optval = x' * ( Q * x ) + v' * x + w;

    else

        %
        % Constant x, affine Q, complex case
        %

        xR = real( x );
        xI = imag( x );
        cvx_optval = xR' * ( real( v ) + real( Q ) * xR ) + xI' * ( imag( v ) + imag( Q ) * xI ) + w;

    end

elseif ~cvx_isaffine( x ),

    error( 'First argument must be affine.' );
    
elseif sQ( 1 ) == 1 && cvx_constant( Q ) ~= 0,
    
    %
    % Constant scalar Q, affine x
    %
    
    cvx_optval = real( Q ) * sum_square_abs( x ) + v' * x + w;
    return
    
else

    %
    % Constant matrix Q, affine x
    %

    cvx_optval = 0;
    while true,
        Q = cvx_constant( Q );
        Q = 0.5 * ( Q + Q' );
        
        %
        % Quick exit for a zero Q
        %

        nnzQ = nnz( Q );
        if nnzQ == 0,
            break
        end
        
        %
        % Remove zero rows and columns from Q. If a diagonal element of Q is
        % zero but there are elements on that row or column that are not,
        % then we know that neither Q nor -Q is PSD.
        %
        
        dQ = diag( Q );
        trQ = sum(dQ);
        if ~all( dQ ),
            tt = dQ ~= 0;
            Q = Q( tt, tt );
            if nnz( Q ) ~= nnzQ,
                success = false;
                break
            end
            dQ = dQ( tt );
            if nnz( v ),
                cvx_optval = v( ~tt, : )' * cvx_subsref( x, ~tt, ':' );
                v = v( tt, : );
            end
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
        Q = 0.5 * ( Q + Q' );
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
            [ R, p ] = chol(full(Q));
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
        cvx_optval = cvx_optval + sg * alpha * sum_square_abs( ( R * x + vbar ) / sqrt(alpha) ) + sum( v' * x ) + wbar;
        break;
        
    end
    
    if ~success,
        error( 'The second argument must be positive or negative semidefinite.' );
    end

end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
