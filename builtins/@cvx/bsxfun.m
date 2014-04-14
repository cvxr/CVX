function z = bsxfun( func, x, y, varargin )

global cvx___
oflag = cvx___.broadcast;
try
    cvx___.broadcast = true;
    z = func( x, y, varargin{:} );
    cvx___.broadcast = oflag;
catch exc
    cvx___.broadcast = oflag;
    if strncmp( exc.identifier, 'CVX:', 4 ), throw( exc ); 
    else rethrow( exc ); end
end


