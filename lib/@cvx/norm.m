function cvx_optval = norm( varargin )
error( nargchk( 1, 3, nargin ) );

%NORM   Matrix or vector norm.
%
%   For CVX, the norm function has been modified to include a single
%   additional case, the 'largest-k' norm, which is introduced below. It
%   has also been adapted to be used in CVX objective functions and
%   expressions. It calls the built-in version of NORM in all other cases
%   to insure identical behavior.
%
%   For vectors...
%     NORM(V,P)           = sum(abs(V).^P)^(1/P)
%     NORM(V)             = norm(V,2).
%     NORM(V,inf)         = max(abs(V)).
%     NORM(V,-inf)        = min(abs(V)).
%   Note that NORM(V,P) for -Inf<=P<1 is NOT, in fact, a valid norm. MATLAB
%   allows them nonetheless, so CVX does as well---but for numeric vectors
%   ONLY. Such "norms" cannot be used within disciplined convex programs.
%
%   For matrices...
%     NORM(X)             = max(svd(X)).
%     NORM(X,2)           = norm(X,2).
%     NORM(X,1)           = max(sum(abs(X))).
%     NORM(X,Inf)         = max(sum(abs(X'))).
%     NORM(X,'fro')       = norm(X(:),2).
%   NORM(X,P) is not implemented for matrices for any other values of P.
%
%   Disciplined convex programming information:
%       NORM is convex, except when P<1, so an error will result if these
%       non-convex "norms" are used within CVX expressions. NORM is
%       nonmonotonic, so its input must be affine.

%
% Quick exit for cvx_constant case
%

x = varargin{1};
if cvx_isconstant( x ),
    cvx_optval = norm( cvx_constant( x ), varargin{2:end} );
    return
elseif ndims( x ) > 2,
    error( 'norm is not defined for N-D arrays.' );
end
[ m, n ] = size( x );
ismat = m ~= 1 & n ~= 1;

%
% Check second argument
%

if nargin < 2,
    p = 2;
else,
    p = varargin{2};
end

%
% Numeric norms
%

if isnumeric( p ),
    if nargin > 2,
        error( 'Too many arguments.' );
    elseif any( size( p ) ~= 1 ) | ~isreal( p ) | p < 1,
        error( 'Second argument must be a real number 1 or larger.' );
    end
    if ~ismat & length( x ) == 1,
        cvx_optval = abs( x );
    else,
        switch p,
            case 1,
                if ismat,
                    cvx_optval = max( sum( abs( x ), 1 ), [], 2 );
                else,
                    cvx_optval = sum( abs( x ) );
                end
            case 2,
                if ismat,
                    cvx_optval = sigma_max( x );
                else,
                    cvx_begin
                        variable z
                        minimize z
                        { x, z } == lorentz( size( x ), [], ~isreal( x ) );
                    cvx_end
                end
            case Inf,
                if ismat,
                    cvx_optval = max( sum( abs( x ), 2 ), [], 1 );
                else,
                    cvx_optval = max( abs( x ) );
                end
            otherwise,
                error( sprintf( 'Sorry, %g-norms not implemented yet.', p ) );
        end
    end
    return
end

%
% String norms
%

if ischar( p ),
    if size( p, 1 ) ~= 1,
        error( 'Second arugment must be a single string or a real number.' );
    elseif strcmpi( p, 'fro' ),
        if nargin > 2,
            error( 'Too many arguments.' );
        else,
            cvx_optval = norm( svec( x, 2 ), 2 );
        end
    else,
        error( [ 'Invalid norm: ', p ] )
    end
    return
end

error( 'The only norms available are 1 <= p <= Inf and "fro".' );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
