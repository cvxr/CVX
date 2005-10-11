function x = subsref( x, S )
error( cvx_subsref_check( nargin, 2, S ) );
switch S(1).type,
    case '.',
        try,
            x = feval( S(1).subs, x );
            S(1) = [];
        catch,
            error( [ 'Invalid field: ', S(1).subs ] );
        end
    case '()',
        error( 'Cannot use subscripts on dual variables.' );
    case '{}',
        error( 'Cell contents reference from a non-cell array object.' );
    otherwise,
        error( 'Invalid subscript structure.' );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
