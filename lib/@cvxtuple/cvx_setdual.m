function x = cvx_setdual( x, y )
x.dual_ = y;
x.value_ = do_setdual( x.value_, y );
function x = do_setdual( x, y )
switch class( x ),
    case 'struct',
        nx = numel( x );
        if nx > 1,
            y(end+1) = struct( 'type', '()', 'subs', {{}} );
            for k = 1 : numel(x),
                y(end).subs{1} = k;
                x(k) = do_setdual( x(k), y );
            end
        else
            y(end+1) = struct( 'type', '{}', 'subs', '' );
            for k = 1 : length(f),
                y(end).subs = {1,k};
                x.(f{k}) = do_setdual( x.(f{k}), y );
            end
        end
    case 'cell',
        y(end+1) = struct( 'type', '{}', 'subs', {{}} );
        for k = 1 : numel(x),
            y(end).subs{1} = k;
            x{k} = do_setdual( x{k}, y );
        end
    case 'cvx',
        x = cvx_setdual( x, y );
    case 'double',
        x = cvx_setdual( cvx( x ), y );
    otherwise,
        cvx_throw( 'Cannot attach a dual variable to an object of type %s.\n', class( x ) );
end


% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
