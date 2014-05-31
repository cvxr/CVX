function cvx_dcp_error( op, mode, x, varargin )

nargs = nargin - 3;
y = [];
switch mode,
case { 'binary', 'mmult' },
    y = varargin{1};
    varargin(1) = [];
    nargs = nargs - 1;
end
strx = {};
maxed = false;
switch mode,
    case 'unary',
        for k = 1 : numel(x),
            if length(strx)==4, maxed = true; break; end
            s = get_arg( x, k );
            if ~any(strcmp(s,strx)), strx{end+1} = s; end %#ok
        end
    case 'binary',
        nx = numel(x); ny = numel(y);
        if nx == 1, sx = get_arg( x, 1 ); end
        if ny == 1, sy = get_arg( y, 1 ); end
        for k = 1 : max(numel(x),numel(y)); 
            if length(strx)==4, maxed = true; break; end
            if nx > 1, sx = get_arg( x, k ); end
            if ny > 1, sy = get_arg( y, k ); end
            if op(1) >= 'a' && op(1) <= 'z',
                s = [ sx, ', ' sy ];
            else
                s = [ sx, ' ', op, ' ', sy ];
            end
            if ~any(strcmp(s,strx)), strx{end+1} = s; end %#ok
        end
    case 'reduce',
        for k = 1 : size(x,2),
            if length(strx)==4, maxed = true; break; end
            s = get_arg( { x(:,k) }, 1 );
            if ~any(strcmp(s,strx)), strx{end+1} = s; end %#ok
        end
    case 'mmult',
        strx{1} = get_arg( { x }, 1 );
    case 'misc',
        if iscell( x ),
            strx = x;
        else
            strx = { x };
        end
end
for k = 1 : nargs,
    strx = strcat( strx, [ ', ', get_arg( varargin, k ) ] );
end
if op(1) >= 'a' && op(1) <= 'z',
    strx = strcat( { [ op, '( ' ] }, strx );
    strx = strcat( strx, ' )' );
elseif ~isequal( mode, 'binary' ),
    strx = strcat( { [ op, ' ' ] }, strx );
end
if length( strx ) == 1,
    plural = '';
    strx = strx{1};
else
    plural = 's';
    if maxed, strx{end+1} = '... and others'; end
    strx = sprintf( '\n        %s', strx{:} );
end
switch op,
    case { '<', '<=', '>', '>=', '==' },
        type = 'constraint';
    case '~=',
        type = 'constraint';
        strx = sprintf( '%s\n   Not-equal (~=) constraints are never allowed in CVX.', strx );
    otherwise,
        type = 'operation';
end
if any(strfind(strx,'nomial')), 
    gtype = 'geometric'; 
else
    gtype = 'convex'; 
end
cvx_throw( 'Disciplined %s programming error:\n   Invalid %s%s: %s', gtype, type, plural, strx );

function strx = get_arg( x, vk )
if ~iscell( x ) || numel( x ) > 1,
    x = x(vk);
else
    x = x{vk};
end
if isempty( x ),
    strx = '[]';
    return
elseif cvx_isconstant( x ),
    x = cvx_constant( x );
end
if isa( x, 'cvx' ) || isnumeric( x ) && numel( x ) > 1,
    strx = [ '{', cvx_class( x, true, true, true ) ];
    if numel( x ) > 1, 
        sx = size( x );
        if length(sx) > 2, tp = 'array';
        elseif all(sx>1), tp = 'matrix';
        else tp = 'vector'; 
        end
        strx = [ strx, ' ', tp ];
    end
    strx = [ strx, '}' ];
elseif ~isnumeric( x ),
    strx = [ '{', class( x ), '}' ];
elseif isreal( x ),
    strx = sprintf( '{%g}', x );
else
    strx = sprintf( '{%g+1j*%g}', real(x), imag(x) );
end
