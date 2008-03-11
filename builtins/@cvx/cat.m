function y = cat( dim, varargin )

%Disciplined convex/geometric programming information for CAT:
%   CAT imposes no convexity restrictions on its arguments.

if ~isnumeric( dim ) | any( size( dim ) ~= 1 ) | dim <= 0 | dim ~= floor( dim ),
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
% Check sizes
%

nz = 0;
msiz = 0;
nzer = 0;
first = 1;
nargs = 0;
msimp = 1;
isr = true;
for k = 1 : nargin - 1,
    x    = cvx( varargin{ k } );
    sx   = x.size_;
    bx   = x.basis_;
    nz   = max( nz, size( bx, 1 ) );
    nzer = nzer + nnz( bx );
    if length( sx ) < dim,
        sx( end + 1 : dim ) = 1;
    end
    if k == 1 | ( first & all( sx ) ),
        lsiz = sx( 1 : dim -1 );
        rsiz = sx( dim + 1 : end );
        first = any( sx == 0 );
    elseif ~all( sx == 0 ),
        if ~isequal( lsiz, sx( 1 : dim - 1 ) ) | ~isequal( rsiz, sx( dim + 1 : end ) ),
            error( sprintf( 'All dimensions but the one being concatenated (%d) must be equal.', dim ) );
        end
    end
    msiz = msiz + sx( dim );
    if all( sx ),
        nargs = nargs + 1;
        varargin{nargs} = x;
        if ~isreal( bx ),
            isr = false;
        end
    end
end

%
% Concatenate
%

sz   = [ lsiz, msiz, rsiz ];
lsiz = prod( lsiz );
rsiz = prod( rsiz );
if nzer == 0,

    yb = [];

elseif nargs == 1,

    y = varargin{1};
    return

elseif lsiz == 1 | rsiz == 1,

    psz = prod( sz );
    issp = cvx_use_sparse( [ nz, psz ], nzer, isr );
    for k = 1 : nargs,
        x = varargin{k};
        x = x.basis_;
        if issp ~= issparse( x ),
            if issp, x = sparse( x ); else x = full( x ); end
        end
        qx = size( x, 1 );
        if qx < nz,
            if issp,
                x = [ x ; sparse( nz - qx, size( x, 2 ) ) ];
            else,
                x = [ x ; zeros( nz - qx, size( x, 2 ) ) ];
            end
        end
        if rsiz > 1,
            x = cvx_reshape( x, [ prod( size( x ) ) / rsiz, rsiz ] );
        end
        varargin{k} = x;
    end
    if rsiz == 1,
        yb = builtin( 'horzcat', varargin{1:nargs} );
    else
        yb = builtin( 'vertcat', varargin{1:nargs} );
        yb = cvx_reshape( yb, [ nz, psz ] );
    end

else

    ndx = 0;
    psz = prod( sz );
    if cvx_use_sparse( [ nz, psz ], nzer, isr ),
        yb = sparse( [], [], [], nz, psz, nzer );
    else
        yb = zeros( nz, psz );
    end
    for k = 1 : nargs,
        x = varargin{ k };
        s = x.size_;
        b = x.basis_;
        s = [ s, ones( 1, dim - length( s ) ) ];
        tsiz = s( dim );
        ltsiz = lsiz * tsiz;
        rndx = ndx + 1 : ndx + ltsiz;
        ndx  = rndx( end );
        if rsiz > 1,
            rndx = rndx( : );
            rndx = rndx( :, ones( 1, rsiz ) ) + nmvec( ones( 1, ltsiz ), : );
        end
        yb( 1 : size( b, 1 ), rndx( : ) ) = b;
    end

end

%
% Create object
%

y = cvx( sz, yb );

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
