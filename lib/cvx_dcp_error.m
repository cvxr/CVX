function estr = cvx_dcp_error( varargin )

x = varargin{1};
estr = '';
if isempty( x ), return; end
op = varargin{end};
y = [];
tt = [];
if nargin == 3,
    if iscell( varargin{2} ),
        x = [ x, varargin{2} ];
    else
        tt = varargin{2};
    end
elseif nargin == 4,
    y = varargin{2};
    tt = varargin{3};
end
if iscell( x ) && size( x, 2 ) == 2,
    y = x(:,2);
    x = x(:,1);
end
if isa( tt, 'cvx' ),
    tt = cvx_classify( tt ) == 19;
end
if ~isempty(tt),
    tt = find(tt);
    if isempty(tt),
        return;
    end
end
strs = {};
for k = 1 : max(numel(x),numel(y)),
    if isempty( tt ), vk = k; else vk = tt(k); end
    strx = get_arg( x, vk );
    if isempty( y ),
        if op(1) >= 'a' && op(1) <= 'z',
            strx = [ op, '( ', strx, ' )' ]; %#ok
        else
            strx = [ op, ' ', strx ]; %#ok
        end
    else
        stry = get_arg( y, vk );
        if op(1) >= 'a' && op(1) <= 'z',
            strx = [ op, '( ', strx, ', ', stry, ' )' ]; %#ok
        else
            strx = [ strx, ' ', op, ' ', stry ]; %#ok
        end
    end
    strs{end+1} = strx; %#ok
end
[strs,ia] = unique(strs);
[ia,ndxs] = sort(ia); %#ok
strs = strs(ndxs);
if length( strs ) == 1,
    plural = '';
    strs = strs{1};
else
    plural = 's';
    strs = sprintf( '\n        %s', strs{:} );
end
switch op,
    case { '<', '<=', '>', '>=', '==' },
        type = 'constraint';
    case '~=',
        type = 'constraint';
        strs = sprintf( '%s\n   Not-equal (~=) constraints are never allowed in CVX.', strs );
    otherwise,
        type = 'operation';
end
try
    error( 'CVX:DCPError', 'Disciplined convex programming error:\n   Invalid %s%s: %s', type, plural, strs );
catch estr
    if nargout == 0,
        rethrow( estr );
    end
end

function strx = get_arg( x, vk )
if ~iscell( x ) || numel( x ) > 1,
    x = x(vk);
else
    x = x{vk};
end
if cvx_isconstant( x ),
    x = cvx_constant( x );
    if numel(x) > 1 && ~nnz( x ~= x(1) ),
        x = x(1);
    end
end
if isa( x, 'cvx' ) || isnumeric( x ) && numel( x ) > 1,
    strx = [ '{', cvx_class( x, true, true, true ), '}' ];
elseif ~isnumeric( x ),
    strx = [ '{', class( x ), '}' ];
elseif isreal( x ),
    strx = sprintf( '%g', x );
else
    strx = sprintf( '%g+1j*%g', real(x), imag(x) );
end
