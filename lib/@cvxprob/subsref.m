function x = subsref( x, S, cheat )

%
% Check the subscript argument
%

if nargin < 3 | isempty( cheat ),
    err = cvx_subsref_check( nargin, 2, S );
    if ~isempty( err ), error( err ); end
end

%
% Process first subscript
%

switch S(1).type,
    case '.',
        if ~ischar( S(1).subs ) | size( S(1).subs, 1 ) > 1,
            error( 'Invalid subscript structure.' );
        end
        switch S(1).subs,
	    case { 'id', 'index' },
            x = feval( S(1).subs, x );
        otherwise,
            try
                global cvx___
                x = subsref( cvx___.problems( index( x ) ), S(1) );
            catch,
                error( [ 'Invalid field: ', S(1).subs ] );
            end
        end
    case '()',
        error( 'Array subscripting not supported for cvx problem objects.' );
    case '{}',
        error( 'Cell subscripting not supported for cvx problem objects.' );
    otherwise,
        error( 'Invalid subscript structure.' );
end

%
% Process remaining subscripts
%

if length( S ) > 1,
    x = subsref( x, S( 2 : end ) );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
