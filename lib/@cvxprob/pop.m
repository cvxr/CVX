function pop( self, clearmode )

global cvx___
cvx_global

try
    p = self.index_;
    prob = cvx___.problems( p );
    if self ~= prob.self, pid = 0;
    else pid = cvx_id( self ); end
catch
    pid = 0;
end

if ~pid,
    p = 1;
    if ~isempty( cvx___.problems ),
        prob = cvx___.problems( 1 );
    else
        prob = [];
    end
end

if ~isempty( prob ) && prob.cleared,
    erase( self ); 
end

s1 = evalin( 'caller', 'who' );
s2 = cellfun( @cvx_id, evalin( 'caller', 'cellfun(@eval,who,''UniformOutput'',false)' ) );
if isequal( clearmode, 'value' ),
    tt = s2 > pid;
    temp = sprintf( '%s,', s1{tt} );
    if ~isempty( temp ),
        temp(end) = [];
        evalin( 'caller', sprintf( '[%s]=cvx_values(%s);', temp, temp ) );
    end
end
if isequal( clearmode, 'clear' ),
    tt = s2 >= pid;
else
    tt = s2 == pid;
end
temp = sprintf( '%s,', s1{tt} );
if ~isempty(s1),
    evalin( 'caller', [ 'clear ', temp(1:end-1) ] );
end

cvx___.problems( p : end ) = [];
cvx___.x = [];
cvx___.y = [];

if p == 1,
    cvx_clearpath( 1 );
    if cvx___.profile,
        profile off;
    end
end

if ~pid,
    error( 'CVX:InternalError', 'Internal CVX data corruption. Please rebuild your model.' );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
