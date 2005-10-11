function x = subsref( x, S, cheat )
error( cvx_verify( x ) );

%
% Check the subscript argument
%

err = cvx_subsref_check( nargin, 2, S );
if ~isempty( err ), error( err ); end

%
% Process first (and possibly second) subscript
%

switch S(1).type,
    case '()',
        if ~isempty( S(1).subs ),
            sz_o = size( x );
            ndxs = reshape( 1 : prod( sz_o ), sz_o );
            try,
                ndxs = ndxs( S(1).subs{:} );
            catch,
                error( lasterr );
            end
            b = cvx_basis( x );
            x = cvx( problem( x ), size( ndxs ), b( ndxs, : ) );
        end
        S(1) = [];
    case '.',
        switch S(1).subs,
            case 'basis',
                sz = size( x );
                x  = cvx_basis( x );
                S(1) = [];
                if ~isempty( S ) & strcmp( S(1).type, '()' ),
                    if isempty( S(1).subs ),
                        S(1) = [];
                    elseif ~isnumeric( S(1).subs{1} ) | length( S(1).subs{1} ) ~= 1,
                        error( 'Basis subscript must be a nonnegative integer.' );
                    else,
                        bndx = S(1).subs{1};
                        if bndx < 0 | bndx ~= floor( bndx ),
                            error( 'Basis subscript must be a nonnegative integer.' );
                        end
                        if bndx >= size( x, 2 ),
                            x = sparse( [], [], [], size( x, 1 ), 1 );
                        else,
                            x = x( :, bndx + 1 );
                        end
                        x = reshape( x, sz );
                        if length( S(1).subs ) == 1,
                            S(1) = [];
                        else,
                            S(1).subs(1) = [];
                        end
                    end
                end
            otherwise,
            	try,
                    x = feval( S(1).subs, x );
            	    S(1) = [];
                catch,
                    error( [ 'Unrecognized field name: ', S(1).subs ] );
                end
        end
    case '{}',
        error( 'Cell contents reference from a non-cell array object.' );
    otherwise,
        error( 'Invalid subscript structure.' );
end

%
% Process remaining subscripts
%

if ~isempty( S ),
    x = subsref( x, S );
end

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
