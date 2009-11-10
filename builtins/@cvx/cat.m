function y = cat( dim, varargin )

%Disciplined convex/geometric programming information for CAT:
%   CAT imposes no convexity restrictions on its arguments.

if ~isnumeric( dim ) || any( size( dim ) ~= 1 ) || dim <= 0 || dim ~= floor( dim ),
    error( 'First argument must be a dimension.' );
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
nzer  = 0;
nargs = 0;
isr   = true;
for k = 1 : nargin - 1,
    x = varargin{k};
    sx = size( x );
    if any( sx ),
        x = cvx( x );
        bx = x.basis_;
        nz = max( nz, size( bx, 1 ) );
        nzer = nzer + nnz( bx );
        sx( end + 1 : dim ) = 1;
        if isempty( sz ),
            sz = sx;
        elseif length( sx ) ~= length( sz ) || nnz( sx - sz ) > 1,
            error( 'All dimensions but the one being concatenated (%d) must be equal.', dim );
        else
            sz( dim ) = sz( dim ) + sx( dim );
        end
        if ~isreal( bx ), 
            isr = false; 
        end
        if all( sx ),
            nargs = nargs + 1;
            varargin{nargs} = x;
        end
    end
end

%
% Simple cases
%

if nargs == 0,
    
    if isempty( sz ), sz( dim ) = nargin; end
    y = cvx( sz, [] );
    return
    
elseif nargs == 1,
    
    y = varargin{1};
    return
    
end

%
% Harder cases
%

msiz = sz( dim );
lsiz = prod( sz( 1 : dim - 1 ) );
rsiz = prod( sz( dim + 1 : end ) );
psz  = lsiz * msiz * rsiz;
if rsiz > 1 && cvx_use_sparse( [ nz, psz ], nzer, isr ),

    cmsiz = 0;
    bi = cell( 1, nargs );
    bj = cell( 1, nargs );
    bv = cell( 1, nargs );
    for k = 1 : nargs,
        x = varargin{k};
        sx = x.size_;
        sx(end+1:dim) = 1;
        [ bi{k}, bj{k}, bv{k} ] = find( x.basis_ );
        bj{k} = bj{k} + lsiz * ( cmsiz + floor((bj{k}-1)/lsiz) * msiz );
        cmsiz = cmsiz + sx(dim);
    end
    bi = builtin( 'vertcat', bi{:} );
    bj = builtin( 'vertcat', bj{:} );
    bv = builtin( 'vertcat', bv{:} );
    yb = sparse( bi, bj, bv, nz, psz );
    
else
    
    for k = 1 : nargs,
        x = varargin{k}.basis_;
        x( end + 1 : nz, end ) = 0;
        if rsiz > 1,
            x = reshape( x, numel(x) / rsiz, 1, rsiz );
        end
        varargin{k} = x;
    end
    yb = builtin( 'cat', 2, varargin{1:nargs} );
    yb = reshape( yb, nz, psz );
    
end
 
%
% Create object
%

y = cvx( sz, yb );

% Copyright 2009 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
