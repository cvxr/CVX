function cvx_clearpath( arg )

% CVX_CLEARPATH   Clears the cvx path.
%
% CVX_CLEARPATH removes the internal cvx directories from Matlab's path. CVX
% does this automatically when a model is completed (i.e., after CVX_END), in
% order to reduce potential naming conflicts with other packages. There is no
% need to call this function during the normal use of CVX.

if nargin == 0,
    cvx_clear( 1 );
end
global cvx___
if ~isempty( cvx___ ) & cvx___.path.active & isempty( cvx___.problems ) & ( nargin == 0 | ~cvx___.path.hold ),
    if ~isempty( cvx___.path.string ),
        opath = matlabpath;
        temp = strfind( opath, cvx___.path.string );
        if ~isempty( temp ),
            opath( temp(1) : temp(1) + length(cvx___.path.string) - 1 ) = [];
        end
        matlabpath( opath );
    end
    cvx___.path.active = false;
end
if nargin == 0,
    cvx___.path.hold = false;
end

% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
