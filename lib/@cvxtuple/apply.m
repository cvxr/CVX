function y = apply( func, x )
y = do_apply( func, x.value_ );

function y = do_apply( func, x )
global cvx___
switch class( x ),
    case 'struct',
        y = cell2struct( do_apply( func, struct2cell( x ) ), fieldnames( x ), 1 );
    case 'cell',
        if cvx___.hcellfun,
            y = cellfun( func, x, 'UniformOutput', false );
        else
            sx = size( x );
            y = cell( sx );
            for k = 1 : prod( sx ),
                y{k} = feval( func, x{k} );
            end
        end
    otherwise,
        y = feval( func, x );
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.


