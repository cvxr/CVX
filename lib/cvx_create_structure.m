function str = cvx_create_structure( sz, varargin )
error( nargchk( 1, Inf, nargin ) );

[ temp, sz ] = cvx_check_dimlist( sz, false );
if ~temp,
    error( 'First argument must be a non-empty dimension list.' );
elseif nargin <= 1,
    nel = prod( sz );
    str = sparse( 1 : nel, 1 : nel, 1, nel, nel );
    return
end

nstrs = nargin - 1;
strs  = cell( 1, nstrs );
sarg  = zeros( 1, nstrs );
rarg  = zeros( 1, nstrs );
for k = 1 : nstrs,
    argt = varargin{k};
    if isstruct( argt ),
        args = eval( 'argt.args', '{}' );
        argt = argt.type;
    else,
        args = {};
    end
    name = [ 'cvx_s_', argt ];
    temp = exist( name );
    if temp ~= 2 & temp ~= 3,
        error( sprintf( 'Undefined matrix structure type: %s', argt ) );
    end
    if iscell( args ),
        str = feval( name, sz( 1 ), sz( 2 ), args{:} );
    else,
        str = feval( name, sz( 1 ), sz( 2 ), args );
    end
    strs{ k } = str;
    sarg( k ) = size( str, 2 );
    rarg( k ) = isreal( str );
end

if length( strs ) == 1,

   strs = strs{1};
   
else,
   
    [ temp, ndxs ] = sort( sarg );
    if any( rarg( 1 : end - 1 ) == rarg( 2 : end ) ),
        temp = min( find( ~rarg( ndxs ) ) );
        if temp < nargin - 1,
            rng  = temp + 1 : nstrs;
            ndx1 = ndxs( rng );
            [ temp2, ndx2 ] = sort( sarg( ndx1 ) .* ( 2 - rarg( ndx1 ) ) );
            ndxs( rng ) = ndx1( ndx2 );
        end
    end

    iscplx = false;
    for k = 1 : nstrs,
        A = strs{ndxs(k)};
        if ~isreal( A ),
            if ~iscplx,
                [ r, c, v ] = find( str );
                r = 2 * r( : ); c = 2 * c( : ); v = v( : );
                str = sparse( [ r - 1 ; r ], [ c - 1 ; c ], [ v, v ], 2 * size( str, 1 ), 2 * size( str, 2 ) );
                iscplx = true;
            end
            [ rr, cr, vr ] = find( real( A ) );
            [ ri, ci, vi ] = find( imag( A ) );
            A = sparse( [ 2 * rr - 1 ; 2 * ri ], [ cr ; ci ], [ vr ; vi ], 2 * size( A, 1 ), size( A, 2 ) );
        elseif iscplx,
            [ r, c, v ] = find( A );
            r = 2 * r( : ); c = 2 * c( : ); v = v( : );
            A = sparse( [ r - 1 ; r ], [ c - 1 ; c ], [ v, v ], 2 * size( A, 1 ), 2 * size( A, 2 ) );
        end
        if k == 1,
            str = A;
        else,
            str = str * cvx_bcompress( str' * A );
        end
        if isempty( str ),
            break;
        end
    end
    
    if iscplx,
        str = str( 1 : 2 : end, : ) + j * str( 2 : 2 : end, : );
    end

end

if length( sz ) > 2,
    str = cvx_replicate_structure( str, sz( 3 : end ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
