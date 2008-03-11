function [ cvx_optval, success ] = quad_form( x, Q, tol )

%QUAD_FORM   Internal cvx version.

%
% Check sizes and types
%

error( nargchk( 2, 3, nargin ) );
if nargin < 3, tol = 4 * eps; end
sx = size( x );
if length( sx ) ~= 2 | all( sx > 1 ),
    error( 'The first argument must be a row or column.' );
else
    sx = prod( sx );
end

sQ = size( Q );
if length( sQ ) ~= 2 | sQ( 1 ) ~= sQ( 2 ),
    error( 'The second argument must be a scalar or a square matrix.' );
elseif sQ( 1 ) ~= sx & sQ( 1 ) ~= 1,
    error( 'Sizes are incompatible.' );
else
    sQ = sQ( 1 );
end

x = vec( x );
success = true;
if cvx_isconstant( x ),

    if isreal( Q ) | isreal( x ),

        %
        % Constant x, affine Q, real case
        %

        x = real( x );
        Q = real( Q );
        cvx_optval = x' * ( Q * x );

    else

        %
        % Constant x, affine Q, complex case
        %

        xR = real( x );
        xI = imag( x );
        cvx_optval = xR' * ( real( Q ) * xR ) + xI' * ( imag( Q ) * xI );

    end

elseif ~cvx_isaffine( x ),

    error( 'First argument must be affine.' );
    
elseif size( Q, 1 ) == 1,
    
    %
    % Constant scalar Q, affine x
    %
    
    cvx_optval = real( Q ) * sum_square_abs( x );
    
else

    %
    % Constant matrix Q, affine x
    %

    cvx_optval = [];
    while true,
        Q = cvx_constant( Q );
        
        %
        % Quick exit for a zero Q
        %

        nnzQ = nnz( Q );
        if nnzQ == 0,
            cvx_optval = 0;
            break
        end
        
        %
        % Remove zero rows and columns from Q. If a diagonal element of Q is
        % zero but there are elements on that row or column that are not,
        % then we know that neither Q nor -Q is PSD.
        %
        
        dQ = diag( Q );
        if ~all( dQ ),
            tt = dQ ~= 0;
            Q = Q( tt, tt );
            if nnz( Q ) ~= nnzQ,
                break
            end
            dQ = dQ( tt );
            x = cvx_subsref( x, tt, ':' );
        end
        
        %
        % Determine the sign of the elements of Q. If they are not all of
        % the same sign, then neither Q nor -Q is PSD.
        %

        dQ = dQ > 0;
        if all( dQ ),
            alpha = +1;
        elseif any( dQ ),
            break
        else
            alpha = -1;
            Q = -Q;
        end
        
        %
        % Now perform a Cholesky factorization. If it succeeds then we
        % know that alpha * Q is PSD.
        %

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
            [ R, p ] = chol( full( Q ) );
            valid = p == 0;
            if ~valid,
                R = [ R, R' \ Q(1:p-1,p:end) ];
            end
        end
        if ~valid,
            valid = false; % normest( Q - R' * R ) < tol * normest( Q );
        end
        
        %
        % If more accuracy is needed, perform a Schur decomposition.
        %
        
        if valid,
            Dmax = max(R(:));
            R = R / Dmax;
            alpha = alpha * Dmax * Dmax;
        else
            [ V, D ] = eig( full( Q ) );
            if cvx_use_sparse( V ),
                V = sparse( V );
            end
            D = diag( D );
            Dmax = max( D );
            Derr = tol * Dmax;
            if min( D ) < - Derr,
                break
            end
            tt = find( D > Derr );
            alpha = alpha * Dmax;
            R = diag(sparse(sqrt(D(tt)/Dmax))) * V( :, tt )';
        end

        cvx_optval = alpha * sum_square_abs( R * x );
        success = true;
        break;
        
    end
    
    if isempty( cvx_optval ),
        if nargout > 1,
            success = false;
        else
            error( 'The second argument must be positive or negative semidefinite.' );
        end
    end

end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
