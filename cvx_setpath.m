function cvx_setpath( arg )

% CVX_SETPATH   Sets the cvx path.
%        CVX_SETPATH adds the internal cvx directories to Matlab's path so
%        that the parser can find the functions they contain. There is no
%        reason to call this function during normal use of cvx; it is done
%        automatically as needed. However, if you are debugging cvx for
%        some reason, this helps to facilitate the process.

% Copyright 2005 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.

% Create the global cvx data structure
global cvx___
if isempty( cvx___ ),
    pstr   = struct( 'string', '', 'formed', false, 'active', false, 'hold', false );
    cvx___ = struct( 'path', pstr, 'problems', [], 'stack', {{}}, 'id', 0, 'has_mex', false, 'pause', false, 'quiet', false );
end

% Set the hold flag
if nargin == 0,
    cvx___.path.hold = true;
end

% Determine the path string to add
if ~cvx___.path.formed,
    opath = matlabpath;
    s = which( 'cvx_begin' );
    if isempty( s ),
        error( 'Cannot set the cvx path unless cvx_begin is in the current path.' );
    end
    if ispc,
        fs = '\';
        ps = ';';
    else
        fs = '/';
        ps = ':';
    end
    s( max( strfind( s, 'cvx_begin' ) ) : end ) = [];
    subs = { 'lib', 'functions', 'sets', 'structures', 'sedumi' };
    if str2num(version('-release')) > 13,
        subs{end+1} = 'keywords';
    end
    npath = '';
    for k = 1 : length( subs ),
        temp = [ s, subs{k} ];
        if exist( temp, 'dir' ),
            temp2 = [ temp, ps ];
            npath = [ npath, temp2 ];
            ndxs = strfind( opath, temp2 );
            if ~isempty( ndxs ),
                opath( ndxs(1) : ndxs(1) + length(temp2) - 1 ) = [];
                matlabpath( opath );
            end
        elseif ~isequal( subs{k}, 'sedumi' ),
            error( [ 'Cannot find the required cvx subdirectory: ', temp ] );
        end
    end
    cvx___.path.string = npath;
end

% Add the string to the path, if it's not already there
if ~cvx___.path.active,
    if ~isempty( cvx___.path.string ),
        matlabpath( [ cvx___.path.string, matlabpath ] );
    end
    if ~cvx___.path.formed,
        cvx___.has_mex = exist( 'cvx_bcompress_mex', 'file' ) == 3;
    end
end

cvx___.path.formed = true;
cvx___.path.active = true;
