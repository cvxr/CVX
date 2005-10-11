function cvx_clearpath( arg )

% CVX_CLEARPATH   Clears the cvx path.
%        CVX_CLEARPATH removes the internal cvx directories to Matlab's 
%        when they are no longer needed, to avoid any possibility that they
%        might shadow the functions from another package. There is no need
%        to call this function during the normal use of cvx; it is done
%        automatically as needed. In fact, if you attempt to use it, and
%        cvx still needs the directories in the path for some reason, it 
%        will simply return silently without doing anything.

if nargin == 0,
    cvx_clear( 1 );
end
global cvx___
if ~isempty( cvx___ ) & cvx___.path.active & isempty( cvx___.stack ) & ( nargin == 0 | ~cvx___.path.hold ),
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

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
