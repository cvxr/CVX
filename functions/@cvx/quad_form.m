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

success = true;
if cvx_isconstant( x ),

    if isreal( Q ) | isreal( x ),

        %
        % Constant x, affine Q, real case
        %

        x = real( x( : ) );
        Q = real( Q );
        cvx_optval = x' * ( Q * x );

    else

        %
        % Constant x, affine Q, complex case
        %

        xR = real( x( : ) );
        xI = imag( x( : ) );
        cvx_optval = xR' * ( real( Q ) * xR ) + xI' * ( imag( Q ) * xI );

    end

elseif ~cvx_isaffine( x ),

    error( 'First argument must be affine.' );
    
elseif size( Q, 1 ) == 1,
    
    %
    % Constant scalar Q, affine x
    %
    
    cvx_optval = real( Q ) * sum_square_abs( x(:) );
    
else,

    %
    % Constant matrix Q, affine x
    %

    x = x( : );
    Q = cvx_constant( Q );
    Q = 0.5 * ( Q + Q' );
    nnzs = find( any( Q ~= 0, 2 ) );
    Q = Q( nnzs, nnzs );
    x = x( nnzs, : );
    
    if all( diag( Q ) <= 0 ),
        alpha = -1;
        Q = -Q;
    else,
        alpha = +1;
    end
    
    if cvx_use_sparse( Q ),
        
        Q = sparse( Q );
        prm = symamd( Q );
        if any( diff( prm ) ~= 1 ),
            x = x( prm, : );
            Q = Q( prm, prm );
        end
        R  = cholinc( Q, 'inf' );
        tt = any( isinf( R ), 2 );
        if any( tt ),
            R( tt, : ) = [];
            valid = normest( Q - R' * R ) < tol * normest( Q );
        else,
            valid = true;
        end
        
    else,
        
        [ R, p ] = chol( full( Q ) );
        valid = p == 0;
        
    end
    
    if ~valid,
        [ V, D ] = eig( full( Q ) );
        D = diag( D );
        Derr = tol * max( D );
        if min( D ) < - Derr,
            if nargout > 1,
                success = false;
                cvx_optval = [];
            else,
                error( 'The second argument must be positive or negative semidefinite.' );
            end
        end
        tt = find( D > Derr );
        R = diag( sqrt( D( tt ) ) ) * V( :, tt )';
    end
    
    cvx_optval = alpha * sum_square_abs( R * x );

end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
