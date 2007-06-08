function s = cvx_pause( flag )

% CVX_PAUSE
%
% CVX_PAUSE(TRUE) instructs CVX to pause and wait for user keypress before and
% after proceeding with the numerical solution of a model. The pauses occur
% within the CVX_END. This is useful for demo purposes.

global cvx___
if isempty( cvx___ ), cvx_setpath( 1 ); end
s = cvx___.pause;
if nargin == 1,
    if length( flag ) ~= 1,
        error( 'Argument must be a numeric or logical scalar.' );
    end
    cvx___.pause = double( flag ) ~= 0;
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
