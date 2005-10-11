function cvx_optval = norms( varargin )
error( nargchk( 1, 4, nargin ) );

%NORMS   Computation of multiple vector norms.
%
%   NORMS( X ) provides a means to compute the norms of multiple vectors
%   packed into a matrix or N-D vector. This is useful for performing
%   max-of-norms or sum-of-norms calculations.
%
%   All of the vector norms, including the false "-inf" norm, supported
%   by NORM() have been implemented in the NORMS() command.
%     NORMS(X,P)           = sum(abs(X).^P).^(1/P)
%     NORMS(X)             = NORMS(X,2).
%     NORMS(X,inf)         = max(abs(X)).
%     NORMS(X,-inf)        = min(abs(X)).
%     NORMS(X,'largest',k) = sum_largest(abs(X),k).
%   If X is a vector, these computations are completely identical to
%   their NORM equivalents. If X is a matrix, a row vector is returned
%   of the norms of each column of X. If X is an N-D matrix, the norms
%   are computed along the first non-singleton dimension.
%
%   NORMS( X, [], DIM ) computes Euclidean norms along the dimension
%   DIM. NORMS( X, P, DIM ) and NORMS( X, 'largest', k, DIM ) compute 
%   their respective norms along the dimension DIM.
%
%   Disciplined convex programming information:
%       NORMS is convex, except when P<1, so an error will result if these
%       non-convex "norms" are used within CVX expressions. NORMS is
%       nonmonotonic, so its input must be affine.

%
% Quick exit for cvx_constant case
%

x = varargin{1};
if cvx_isconstant( x ),
    cvx_optval = norms( cvx_constant( x ), varargin{2:end} );
    return
end
sx = size( x );

%
% Check second argument
%

if nargin < 2 | isempty( varargin{2} ),
    p = 2; 
else,
    p = varargin{2};
end

%
% Numeric norms
%

if isnumeric( p ),
    if nargin > 3,
        error( 'Too many arguments.' );
    elseif any( size( p ) ~= 1 ) | ~isreal( p ) | p < 1,
        error( 'Second argument must be a real number 1 or larger.' );
    end
    if nargin < 3 | isempty( varargin{3} ),
        dim = cvx_default_dimension( sx );
    elseif ~cvx_check_dimension( varargin{3}, false ),
        error( 'Third argument must be a valid dimension.' );
    else,
        dim = varargin{3};
    end
    if length( sx ) < dim | sx( dim ) == 1,
        cvx_optval = abs( x );
    else,
        switch p,
            case 1,
	        cvx_optval = sum( abs( x ), dim );
            case 2,
                sz = sx;
                sz( dim ) = 1;
                cvx_begin
                    variable z
                    minimize z
                    { x, z } == lorentz( sx, dim, ~isreal( x ) );
                cvx_end
            case Inf,
                cvx_optval = max( abs( x ), [], dim );
            otherwise,
                error( sprintf( sprintf( 'Sorry, %g-norms not implemented yet.', p ) ) );
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
    elseif strcmpi( p, 'largest' ),
        if nargin < 3,
            error( 'Not enough arguments.' );
        end
        k = varargin{3};
        if ~isnumeric( k ) | ~isreal( k ),
            error( 'Third argument must be a scalar.' );
        end
        if nargin < 4 | isempty( varargin{4} ),
            dim = cvx_default_dimension( sx );
        elseif ~cvx_check_dimension( varargin{4}, false ),
            error( 'Fourth argument must be a valid dimension.' );
        else,
            dim = varargin{4};
        end
        cvx_optval = sum_largest( abs( x ), k, dim );
    else,
        error( [ 'Invalid vector norm: ', p ] )
    end
    return
end

error( 'The only norms available are 1 <= p <= Inf and "largest".' );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
