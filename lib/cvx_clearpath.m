function cvx_clearpath( arg )

%CVX_CLEARPATH   Clears the cvx path.
%   CVX_CLEARPATH removes the internal cvx directories from Matlab's path. CVX
%   does this automatically when a model is completed (i.e., after CVX_END), in
%   order to reduce potential naming conflicts with other packages. There is no
%   need to call this function during the normal use of CVX.

cvx_global
if nargin == 0,
    cvx___.path.hold = false;
end
if cvx___.path.hold,
    if isempty( cvx___.problems ),
        nsolv = cvx___.solver;
    else
        nsolv = cvx___.problems(end).solver;
    end
    cvx_setspath( nsolv );
elseif cvx___.path.active && isempty( cvx___.problems ),
    cvx_setspath( '' );
    if ~isempty( cvx___.path.string ),
        cpath = matlabpath;
        temp = strfind( cpath, cvx___.path.string );
        if ~isempty(temp),
            cpath(temp(1):temp(1)+length(cvx___.path.string)-1) = [];
            s = warning('off');
            matlabpath(cpath);
            warning(s);
        end
    end
    cvx___.path.active = false;
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
