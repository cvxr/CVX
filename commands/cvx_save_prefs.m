function cvx_save_prefs( mode )

%CVX_SAVE_PREFS   Saves current CVX settings for future MATLAB sessions.
%   CVX_SAVE_PREFS saves the the current global CVX settings to a special
%   prefences file (stored in the "prefdir" directory). This enables CVX to
%   remember your preferred settings (solver, precision, etc.) between
%   MATLAB sessions.

global cvx___
if nargin < 1, mode = 0; end
if ~isfield( cvx___, 'pfile' ),
    cvx_throw( 'CVX is not currently loaded; there are no preferences to save.' );
end
if ~mode,
    fprintf( 'Saving CVX preferences...' );
end
outp.path      = cvx___.path;
outp.license   = cvx___.license;
outp.solvers   = cvx___.solvers;
outp.broadcast = cvx___.broadcast;
outp.expert    = cvx___.expert;
outp.precision = cvx___.precision;
outp.precflag  = cvx___.precflag;
outp.quiet     = cvx___.quiet;
outp.profile   = cvx___.profile;
outp.solvers.map.default = cvx___.solvers.selected;
[ outp.solvers.list.solve, outp.solvers.list.eargs ] = deal( {} );
outp.solvers.active = 0;
if mode == 2,
    pfile = [ cvx___.where, cvx___.fs, 'cvx_prefs.mat' ];
else
    pfile = cvx___.pfile;
    outg = cvx___.gprefs;
    if ~isempty( outg ),
        savep = outp;
        try
            for ff = fieldnames( outp )',
                f = ff{1};
                if isequal( outp.(f), outg.(f) ),
                    outp.(f) = [];
                end
            end
        catch
            outp = savep; %#ok
        end
    end
end
try
    save(pfile,'-struct','outp');
    if ~mode,
       fprintf('done.\n');
    end
catch exc
    if ~mode,
        fprintf( 'FAILED.\n' );
    end
    rethrow( exc );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
