function cvx_setspath( nsolv )
global cvx___

%CVX_SETPATH   Sets the cvx solver path.
%   CVX_SETPATH adds the internal cvx solver directories to Matlab's path
%   so that the CVX system can find the functions that they contain. There 
%   is no reason to call this function during normal use of CVX; it is done
%   automatically as needed. However, if you are debugging CVX, calling
%   this function can help to insure that breakpoints stay valid.

cvx_global
nsolv = lower(nsolv);
osolv = cvx___.path.sactive;
if ~strcmp( nsolv, osolv ),
    needupd = false;
    cpath = matlabpath;
    if ~isempty( osolv ),
        tstr = cvx___.path.solvers.(osolv);
        if ~isempty( tstr ),
            temp = strfind( cpath, tstr );
            if ~isempty(temp),
                cpath(temp(1):temp(1)+length(tstr)-1) = [];
                needupd = true;
            end
        end
    end
    if ~isempty( nsolv ),
        if ~isempty( nsolv ),
            tstr = cvx___.path.solvers.(nsolv);
            if ~isempty( tstr ),
                cpath = [ tstr, cpath ];
                needupd = true;
            end
        end
    end
    if needupd,
        s = warning('off');
        matlabpath(cpath);
        warning(s);
    end
    cvx___.path.sactive = nsolv;
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
