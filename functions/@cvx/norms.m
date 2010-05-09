function y = norms( x, p, dim )

%NORMS   Internal cvx version.

%
% Check second argument
%

error( nargchk( 1, 3, nargin ) );
if nargin < 2 || isempty( p ),
    p = 2;
elseif ~isnumeric( p ) || numel( p ) ~= 1 || ~isreal( p ),
    error( 'Second argument must be a real number.' );
elseif p < 1 || isnan( p ),
    error( 'Second argument must be between 1 and +Inf, inclusive.' );
end

%
% Check third argument
%

sx = size( x );
if nargin < 3 || isempty( dim ),
    dim = cvx_default_dimension( sx );
elseif ~cvx_check_dimension( dim, false ),
    error( 'Third argument must be a valid dimension.' );
end
if isempty( x ) || length( sx ) < dim || sx( dim ) == 1,
    p = 1;
end

%
% Type check
%

persistent remap1 remap2
if isempty( remap2 ),
    remap1 = cvx_remap( 'constant', 'log-convex' );
    remap2 = cvx_remap( 'affine', 'log-convex' );
end
xc = reshape( cvx_classify( x ), sx );
if ~all( remap2( xc( : ) ) ),
    error( 'Disciplined convex programming error:\n   Invalid computation: norms( {%s}, ... )', cvx_class( x, true, true ) );
end

%
% Compute norms
%

switch p,
    case 1,    
        y = sum( abs(x), dim );
    case Inf,  
        y = max( abs(x), [], dim );
    otherwise,
        nd = length( sx );
        nx = sx( dim );
        sy = sx;
        sy( dim ) = 1;
        tt = all( remap1( xc ), dim );
        if all( tt( : ) ),
            y = sum( abs(x) .^ p, dim ) .^ (1/p);
        elseif any( tt( : ) ),
            if dim ~= 1,
                perm = [ dim, 1:dim-1, dim+1:length(sx) ];
                x  = permute( x,  perm );
                tt = permute( tt, perm );
                sy = sy( perm );
            else
                perm = [];
            end
            nv = prod( sy );
            x  = reshape( x, nx, nv );
            y  = cvx( [ 1, nv ], [] );
            xt = cvx_subsref( x, ':', tt );
            y  = cvx_subsasgn( y, tt, sum( (abs(xt)).^p, 1 ) .^ (1/p) );
            tt = ~tt;
            xt = cvx_subsref( x, ':', tt );
            y  = cvx_subsasgn( y, tt, norms( xt, p, 1 ) );
            y  = reshape( y, sy );
            if ~isempty( perm ),
                y = ipermute( y, perm );
            end
        elseif p == 2,
            cvx_begin
                epigraph variable y( sy )
                { cvx_accept_convex(x), y } == lorentz( sx, dim, ~isreal( x ) ); %#ok
            cvx_end
        else
            sw = sx;
            sw(nd+1) = 2;
            cvx_begin
                variable z( sx )
                epigraph variable y( sy )
                if isreal(x), cmode = 'abs'; else cmode = 'cabs'; end
                { cat( nd+1, z, cvx_expand_dim( y, dim, nx ) ), cvx_accept_convex(x) } ...
                    == geo_mean_cone( sw, nd+1, [1/p,1-1/p], cmode ); %#ok
                sum( z, dim ) == y;
            cvx_end
        end
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
