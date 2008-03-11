function [ map, pn, pd, cplx ] = cvx_geomean_map( w )

%CVX_GEOMEAN_MAP    Helper function for geomean->SOCP conversion.
%   This is an internal CVX function involved in the conversion of geomean()
%   function calls to SOCP-solvable form.
%
%   This is an internal CVX function, and as such no checking is performed to
%   insure that the arguments are valid.

w = reshape( w, 1, numel(w) );
base = length( w );
if base == 1,
    p = w;
    [ pn, pd ] = rat( w );
    if pn < 0,
        pnd  = pd - pn;
        pnd2 = pow2(ceil(log2(pnd)));
        pd2  = round( pnd2 / ( 1 - w ) );
        pn2  = pnd2 - pd2;
        if abs( pn / pd - w ) >= abs( pn2 / pd2 - w ),
            pn = pn2; pd = pd2;
        end
        w = [ - pn, pd ];
    elseif pn < pd,
        pd2 = pow2(ceil(log2(pd)));
        pn2 = round(pd2*w);
        if abs( pn / pd - w ) >= abs( pn2 / pd2 - w ),
            pn = pn2; pd = pd2;
        end
        w = [ pn, pd - pn ];
    else
        pn2 = pow2(ceil(log2(pn)));
        pd2 = round(pn2/w);
        if abs( pn / pd - w ) >= abs( pn2 / pd2 - w ),
            pn = pn2; pd = pd2;
        end
        w = [ pd, pn - pd ];
    end
end
wsum = sum( w );
for k = unique( factor( wsum ) ),
    while all( rem( w, k ) == 0 ),
        w = w / k;
        wsum = wsum / k;
    end
end
nx   = length( w );
nq   = nextpow2(wsum);
nw   = pow2(nq);
ndxs = 1 : nx;
if nw > wsum,
    ndxs(end+1) = 0;
    w(end+1) = nw - wsum;
end
%
% Look for candidates with common factors that can be
% extracted to save LMI blocks: e.g.,
%    ( x1^3 x2^2 x2^3 )^(1/8) =
%       = ((sqrt(x1 x3))^(1/6) x2^2)^(1/8)
%
n3 = nx;
map = [];
[ ff, ee ] = log2( w );
ndx1 = find( ff ~= 0.5 & ff ~= 0 );
while length( ndx1 ) > 1,
    % Build cross matrix
    nv = length( ndx1 );
    ww = w( ndx1 );
    ww = ones( nv, 1 ) * ww;
    [ wi, wj, ww ] = find( bitand( tril(ww,-1)', triu(ww,1) ) );
    [ ff, ee ] = log2( ww );
    ndx2 = find( ff ~= 0.5 );
    if isempty( ndx2 ), break; end
    % Greedy: select the largest overlap
    ee = ee(ndx2);
    [ wc_t, wnm ] = max( sum(dec2bin(ww(ndx2))-'0',2)+(1-ee/(max(ee)+1)) );
    % Construct a 2-element geomean
    wi_t     = ndx1(wi(ndx2(wnm)));
    wj_t     = ndx1(wj(ndx2(wnm)));
    n3       = n3 + 1;
    ndxs     = [ ndxs, n3 ];
    map      = [ map, [ ndxs(wi_t) ; ndxs(wj_t) ; n3 ] ];
    % Update the weights
    wt       = bitand( w(wi_t), w(wj_t) );
    w(end+1) = 2 * wt;
    wt       = bitcmp( wt, nq );
    w(wi_t)  = bitand( w(wi_t), wt );
    w(wj_t)  = bitand( w(wj_t), wt );
    % Update the count
    ndx1 = [ ndx1, length(ndxs) ];
    [ ff, ee ] = log2( w(ndx1) );
    ndx1 = ndx1( ff ~= 0.5 & ff ~= 0 );
end
%
% Now do standard left-to-right combining
%    x1^3 x2^2 x3^3 = x1 x1 x1 x2 x2 x3 x3 x3
%       = (x1 x1)(x1 x2)(x2 x3)(x3 x3)
%
for k = 1 : nq,
    tt  = rem( w, 2 ) ~= 0;
    w   = floor( 0.5 * w );
    n12 = ndxs( tt );
    ntt = 0.5 * length( n12 );
    if ntt >= 1,
        n3    = n3(end) + 1 : n3(end) + ntt;
        map   = [ map, [ reshape( n12, 2, ntt ) ; n3 ] ];
        ndxs  = [ ndxs, n3 ];
        w     = [ w, ones( 1, ntt ) ];
    end
end
if ~isempty( map ),
    map(map==0) = map(end);
else
    map = zeros(3,0);
end
cplx = size( map, 2 );
if base == 1 & ( base * cplx > cvx_power_warning ),
    warning( sprintf( [ ...
        'Using x^%g = x^(%g/%g) in a CVX model is about %gx more\n',...
        'complex than using x^2. See the user guide under the topic "rational\n',...
        'powers" for more details---including how to disable this warning.' ], ...
        p, pn, pd, cplx ) );
end

% Copyright 2008 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.



