function cvx_setpath( arg )
global cvx___

%CVX_SETPATH   Sets the cvx path.
%   CVX_SETPATH adds the internal cvx directories to Matlab's path so that the
%   CVX system can find the functions that they contain. There is no reason to 
%   call this function during normal use of CVX; it is done automatically as
%   needed. However, if you are debugging CVX, calling this function can help to
%   insure that breakpoints stay valid.

% Set the hold flag
cvx_global
if ~cvx___.path.active,
    if ~isempty( cvx___.path.string ),
        s = warning('off');
        matlabpath([cvx___.path.string,matlabpath]);
        warning(s);
    end
    cvx___.path.active = true;
end
if nargin == 0 || cvx___.path.hold,
    cvx___.path.hold = true;
    if isempty(cvx___.problems),
        nsolv = cvx___.solver;
    else
        nsolv = cvx___.problems(end).solver;
    end
else
    nsolv = '';
end
cvx_setspath(nsolv);

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
