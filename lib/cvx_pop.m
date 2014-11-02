function cvx_pop( p, pid )

global cvx___
try
    np = length( cvx___.problems );
catch
    np = 0;
end
if nargin == 0,
    p = np;
    if p == 0,
        if np == 0, return; end
    elseif p > np,
        return
    else
        pstr = cvx___.problems( p );
        pid = pstr.id;
    end
elseif p > np
    pid = 0;
elseif p == 0
    pid = -1;
else
    pstr = cvx___.problems( p );
    if pstr.id ~= pid, pid = 0; end
end
s1 = evalin( 'caller', 'who' );
s2 = evalin( 'caller', 'cellfun(@eval,who,''UniformOutput'',false)' );
s2 = cellfun( @cvx_id, s2 );
do_timer = false;
if pid <= 0,
    do_erase = ~isempty( cvx___.problems );
    p = +do_erase;
    tt = s2 >= 0;
elseif ~pstr.finished
    tt = s2 >= pid;
    do_erase = true;
elseif isnan( pstr.direction )
    tt = s2 == pid;
    do_erase = false;
elseif pstr.checkpoint(4) >= pstr.checkpoint(1),
    tt = s2 > pid;
    temp = sprintf( '%s,', s1{tt} );
    if ~isempty( temp ),
        temp = temp(1:end-1);
        evalin( 'caller', sprintf( '[%s]=cvx_values(%s);', temp, temp ) );
    end
    do_timer = length(pstr.tictime)==2;
    tt = s2 == pid;
    do_erase = true;
else
    tt = s2 >= pid;
    do_erase = false;
    if ~isempty( pstr.objective ),
        tt = tt & ( s2 ~= cvx_id( pstr.objective ) );
    end
end
s1 = s1(tt);
if ~isempty(s1)
    temp = sprintf( ' %s', s1{:} );
    evalin( 'caller', [ 'clear', temp ] );
end
if do_erase
    cvx_erase( p );
end
if p > 1
    cvx___.problems( p : end ) = [];
else
    cvx___.problems = [];
    cvx_clearpath( 1 );
end
if pid == 0
    cvx_throw( 'Internal CVX data corruption. Please rebuild your model.' );
end
if do_timer
    tfin = tic;
    try
        tval = tfin - pstr.tictime(1);
        cvx___.timers(2) = cvx___.timers(2) + tval;
        cvx___.timers(3) = cvx___.timers(3) + ( tfin - pstr.tictime(2) );
        for q = 1 : p - 1,
            cvx___.problems(q).tictime = cvx___.problems(q).tictime + tval;
        end
    catch
        tval = double(tfin) - double(pstr.tictime(1));
        cvx___.timers(2) = double(cvx___.timers(2)) + double(tval);
        cvx___.timers(3) = double(cvx___.timers(3)) + ( double(tfin) - double(pstr.tictime(2)) );
        for q = 1 : p - 1,
            cvx___.problems(q).tictime = double(cvx___.problems(q).tictime) + tval;
        end
    end
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
