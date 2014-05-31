function cvx_throw( exc, varargin )
if ischar( exc ),
    if strncmp( exc, 'CVX:', 4 )
        ndxs = find( exc == ':', 2, 'first' );
        e2.identifier = exc(1:ndxs(end)-1);
        exc = exc(ndxs(end)+1:end);
    elseif strncmp( exc, 'Disciplined ', 12 ),
        e2.identifier = 'CVX:DCPError';
    else
        e2.identifier = 'CVX:Error';
    end
    if nargin == 1,
        e2.message = sprintf( exc, '' );
    else
        e2.message = sprintf( exc, varargin{:} );
    end
elseif ~strncmp( exc.identifier, 'CVX:', 4 ),
    rethrow( exc );
else
    e2 = struct( 'identifier', exc.identifier, 'message', exc.message );
end
e2.stack = dbstack;
if ischar( exc ),
    ns = find(~strncmp({e2.stack.name},'cvx_',4),1,'first');
else
    ns = min(2,length(e2));
end
e2.stack = e2.stack(ns:end);
error( e2 );
