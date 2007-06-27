function y = cvx_id( x )
y = apply( @cvx_id, x );
switch class( y ),
    case 'struct',
        y = struct2cell( y );
        y = max( [ -Inf, y{:} ] );
    case 'cell',
        y = max( [ -Inf, y{:} ] );
end

