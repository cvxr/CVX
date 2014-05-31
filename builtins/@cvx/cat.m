function y = cat( dim, varargin )

%Disciplined convex/geometric programming information for CAT:
%   CAT imposes no convexity restrictions on its arguments.

if ~isnumeric( dim ) || any( size( dim ) ~= 1 ) || dim <= 0 || dim ~= floor( dim ),
    cvx_throw( 'First argument must be a dimension.' );
end

%
% Quick exit
%

if nargin == 2,
    y = varargin{1};
    return
end

%
% Determine the final size and check consistency
%

sz    = [];
nz    = 1;
nargs = 0;
for k = 1 : nargin - 1,
    x = varargin{k};
    sx = size( x );
    if any( sx ),
        bx = cvx_basis( x );
        varargin{k} = bx;
        nz = max( nz, size( bx, 1 ) );
        sx( end + 1 : dim ) = 1;
        if isempty( sz ),
            sz = sx;
        elseif length( sx ) ~= length( sz ) || nnz( sx - sz ) > 1,
            cvx_throw( 'All dimensions but the one being concatenated (%d) must be equal.', dim );
        else
            sz( dim ) = sz( dim ) + sx( dim ); %#ok
        end
        if all( sx ),
            nargs = nargs + 1;
            varargin{nargs} = bx;
        end
    end
end

%
% Simple cases
%

if nargs == 0,
    
    if isempty( sz ), 
        sz( dim ) = nargin; 
    end
    y = cvx( sz, [] );
    return
    
elseif nargs == 1,
    
    y = cvx( sz, varargin{1} );
    return
    
end

%
% Harder cases
%

msiz = sz( dim );
lsiz = prod( sz( 1 : dim - 1 ) );
rsiz = prod( sz( dim + 1 : end ) );
psz  = lsiz * msiz * rsiz;
for k = 1 : nargs,
    x = varargin{k};
    x( end + 1 : nz, end ) = 0;
    if rsiz > 1,
        x = reshape( x, numel(x) / rsiz, rsiz );
    end
    varargin{k} = x;
end
if rsiz > 1,
    yb = builtin( 'cat', 1, varargin{1:nargs} );
else
    yb = builtin( 'cat', 2, varargin{1:nargs} );
end
yb = reshape( yb, nz, psz );
 
%
% Create object
%

y = cvx( sz, yb );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
