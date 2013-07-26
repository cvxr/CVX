function cvx_save_prefs( in_setup )

%CVX_SAVE_PREFS   Saves current CVX settings for future MATLAB sessions.
%   CVX_SAVE_PREFS saves the the current global CVX settings to a special
%   prefences file (stored in the "prefdir" directory). This enables CVX to
%   remember your preferred settings (solver, precision, etc.) between
%   MATLAB sessions.

if nargin < 1, in_setup = false; end
global cvx___
if isempty( cvx___ ), return; end
cvx___.solvers.map.default = cvx___.solvers.selected;
osolv = cvx___.solvers;
try
    [ cvx___.solvers.list.check, cvx___.solvers.list.solve, cvx___.solvers.list.eargs ] = deal( {} );
    cvx___.solvers.active = 0;
    if strncmp( computer, 'PC', 2 ), fs = '\'; fsre = '\\'; else fs = '/'; fsre = '/'; end
    pfile = [ regexprep( prefdir(1), [ fsre, 'R\d\d\d\d\w$' ], '' ), fs, 'cvx_prefs.mat' ];
    save(pfile,'-struct',...
        'cvx___','expert','precision','precflag',...
        'rat_growth','path','license','solvers');
    cvx___.solvers = osolv;
    pfile2 = [ prefdir, fs, 'cvx_prefs.mat' ];
    if ~strcmp( pfile, pfile2 ) && exist( pfile, 'file' ) && exist( pfile2, 'file' ),
        try delete( pfile2 ); catch end %#ok
    end
catch errmsg
    cvx___.solvers = osolv;
    if in_setup,
        rethrow( errmsg );
    else
        errmsg = cvx_error( errmsg, 67, false, '    ' );
        error( 'CVX:BadPrefsSave', 'CVX encountered the following error attempting to save your preferences:\n%sYour preferences will revert once you exit this MATLAB session.', errmsg  );
    end
end

% Copyright 2005-2013 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
