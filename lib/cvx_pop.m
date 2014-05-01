function cvx_pop( p, do_clear, do_warn )

global cvx___

if isa( p, 'cvxprob' ),
    [ p, pstr ] = verify( p, false );
elseif p > length( cvx___.problems ),
    p = 0;
else
    pstr = cvx___.problems( p );
end
    
if p == 0,
    pid = 0;
    values = false;
    do_clear = true;
else
    if nargin >= 3 && do_warn && ( ~isempty(pstr.objective) || ~isempty(pstr.variables) || ~isempty(pstr.duals) || nnz(pstr.t_variable) > 1 ),
        warning( 'CVX:Empty', 'A non-empty cvx problem already exists in this scope.\n   It is being overwritten.', 1 ); %#ok
    end
    if nargin < 2,
        do_clear = false;
    end
    pid = cvx_id( pstr.self );
    values = pstr.cleared && ~do_clear;
    if do_clear || ~values && nnz( pstr.t_variable ) > 1,
        erase( pstr.self ); 
    end
end

s1 = evalin( 'caller', 'who' );
s2 = evalin( 'caller', 'cellfun(@eval,who,''UniformOutput'',false)' );
s2 = cellfun( @cvx_id, s2 );
if values,
    temp = sprintf( '%s,', s1{s2>pid} );
    if ~isempty( temp ),
        temp = temp(1:end-1);
        evalin( 'caller', sprintf( '[%s]=cvx_values(%s);', temp, temp ) );
    end
    tt = s2 == pid;
elseif p == 0 || nargin >= 2 && do_clear,
    tt = s2 >= pid;
else
    tt = s2 == pid;
end
temp = sprintf( ' %s', s1{tt} );
if ~isempty(s1),
    evalin( 'caller', [ 'clear', temp ] );
end

if p > 1,
    cvx___.problems( p : end ) = [];
else
    cvx___.problems = [];
    cvx_clearpath( 1 );
end
cvx___.x = [];
cvx___.y = [];

if ~pid,
    error( 'CVX:InternalError', 'Internal CVX data corruption. Please rebuild your model.' );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
