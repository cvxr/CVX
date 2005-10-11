function y = cat( dim, varargin )
error( cvx_verify( varargin{:} ) );
[ prob, varargin{ : } ] = cvx_operate( [], varargin{ : } );

%
% Check dimension
%

if ~isnumeric( dim ) | any( size( dim ) ~= 1 ) | dim <= 0 | dim ~= floor( dim ),
    error( 'First argument must be a dimension.' );
end

%
% Check sizes
%

nz = 0;
msiz = 0;
nzer = 0;
first = 1;
for k = 1 : nargin - 1,
    x  = varargin{ k };
    sx = size( x );
    bx = cvx_basis( x );
    if length( sx ) < dim, 
        sx( end + 1 : dim ) = 1; 
    end
    nzer = nzer + nnz( bx );
    nz = max( nz, size( bx, 2 ) );
    sx2 = sx; 
    sx2( dim ) = [];
    if k == 1 | ( first & all( sx > 0 ) ),
        lsiz = sx( 1 : dim -1 );
        rsiz = sx( dim + 1 : end );
        first = any( sx == 0 );
    elseif ~all( sx == 0 ),
        if ~isequal( lsiz, sx( 1 : dim - 1 ) ) | ~isequal( rsiz, sx( dim + 1 : end ) ),
            error( sprintf( 'All dimensions but the one being concatenated (%d) must be equal.', dim ) );
        end
    end
    msiz = msiz + sx( dim );
end

%
% Concatenate
%

sz = [ lsiz, msiz, rsiz ];
lsiz = prod( lsiz );
rsiz = prod( rsiz );
if rsiz > 1,
    nmvec = lsiz * msiz;
    nmvec = 0 : nmvec : nmvec * ( rsiz - 1 );
end
ndx = 0;
y = sparse( [], [], [], prod( sz ), nz, nzer );
if nzer ~= 0,
    for k = 1 : nargin - 1,
        x = varargin{ k };
        s = size( x );
        if all( s ~= 0 ),
            b = cvx_basis( x );
            s = [ s, ones( 1, dim - length( s ) ) ];
            tsiz = s( dim );
            ltsiz = lsiz * tsiz;
            rndx = ndx + 1 : ndx + ltsiz;
            ndx  = rndx( end );
            if rsiz > 1,
                rndx = rndx( : );
                rndx = rndx( :, ones( 1, rsiz ) ) + nmvec( ones( 1, ltsiz ), : );
            end
            y( rndx( : ), 1 : size( b, 2 ) ) = b;
        end
    end
end

%
% Create object
%

y = cvx( prob, sz, y );

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
