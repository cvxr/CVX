function cvx_global

%CVX_GLOBAL   Create CVX's global internal data structure, if needed.
%   CVX_GLOBAL creates a hidden structure CVX needs to do its work. It is
%   harmless for the user to call it, but it is also useless to do so.

global cvx___ 
if isfield( cvx___, 'problems' ), return; end
tstart = tic;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the global data structure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cvx_version(1);

commands = { 'cvx_begin', 'cvx_clear', 'cvx_end', 'cvx_expert', ...
    'cvx_pause', 'cvx_power_warning', 'cvx_precision', 'cvx_profile', ...
    'cvx_quiet', 'cvx_save_prefs', 'cvx_tic', 'cvx_toc' };
c_type = cell(1,length(commands));
[ c_type{:} ] = deal('C');
keywords = { 'in', 'dual', 'epigraph', 'expression', 'expressions', ...
    'hypograph', 'maximize', 'maximise', 'minimize', 'minimise', ...
    'subject', 'variable', 'variables', ...
    'epigraph_', 'hypograph_' };
k_type = cell(1,length(keywords)); 
[ k_type{:} ] = deal('K');
structures = { 'banded', 'binary', 'complex', 'diagonal', 'hankel', ...
    'hermitian', 'integer', 'lower_bidiagonal', 'lower_hessenberg', ...
    'lower_triangular', 'nonnegative', 'scaled_identity', ...
    'skew_symmetric', 'semidefinite', 'sparse', 'symmetric', ...
    'toeplitz', 'tridiagonal', 'upper_bidiagonal', 'upper_hankel', ...
    'upper_hessenberg', 'upper_triangular', 'nonnegative_', ...
    'complex_if', 'hermitian_if', 'semicontinuous', 'semiinteger' };
s_type = cell(1,length(structures)); 
[ s_type{:} ] = deal('S');
reserved = cell2struct( [ c_type, k_type, s_type ], ...
    [ commands, keywords, structures ], 2 );

cvx___.reswords    = reserved;
cvx___.problems    = [];
cvx___.logarithm   = sparse(1,1);
cvx___.exponential = sparse(1,1);
cvx___.classes     = int8(3);
cvx___.cones       = struct( 'type', {}, 'indices', {}, 'slacks', {} );
cvx___.equalities  = cell(1,0);
cvx___.inequality  = false(1,0);
cvx___.n_equality  = 0;
cvx___.x           = zeros( 0, 1 );
cvx___.y           = zeros( 0, 1 );
cvx___.warmstart   = {false};
cvx___.id          = 0;
cvx___.obj         = cvx( struct( 'size_', [0,1], 'basis_', sparse(1,0), 'dual_', '', 'id_', 0 ) );
cvx___.pobj        = cvxprob( struct( 'index_', 0, 'id_', 0 ) );
cvx___.timers      = uint64([tstart,0,0,0]);
cvx___.increment   = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run each shim to connect/reconnect the solvers %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cur_d = pwd;
osolvers = cvx___.solvers.list;
nsolv = length(osolvers);
nrej = 0;
for k = 1 : length(osolvers),
    tsolv = osolvers(k);
    try
        cd(tsolv.spath);
        tsolv.warning = '';
        tsolv = feval( tsolv.sname, tsolv );
    catch errmsg
        if isempty( tsolv.name ),
            tsolv.name = [ tsolv.spath, tsolv.sname ];
        end
        tsolv.error = errmsg;
    end
    if ~isempty(tsolv.error) || isempty(tsolv.solve),
        nrej = nrej + 1;
    end
    try
        if k == 1, solvers = tsolv;
        else solvers = [ solvers, tsolv ]; end %#ok
    catch
        for ff = fieldnames(tsolv)',
            solvers(k).(ff{1}) = tsolv.(ff{1});
        end
    end
end
clear osolvers
cvx___.solvers.list = solvers;
cd( cur_d );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If any solvers have errors, force the user to re-run cvx_setup. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~nrej, return; end
reject = {};
ndefault = 0;
for k = 1 : nsolv,
    tsolv = solvers(k);
    terr = tsolv.error;
    if isempty( terr ) 
        if isempty( tsolv.solve ),
            temp = { [ tsolv.name, ': not found' ] };
        else
            temp = {};
            if ndefault == 0, ndefault = k; end
        end
    elseif ~ischar( terr )
        temp = cvx_error( terr, '', [ tsolv.name, ': UNEXPECTED ERROR ' ] );
    elseif any( regexp( terr, '\n' ) )
        temp = cvx_error( terr, '', [ tsolv.name, ': ' ] );
    elseif ~isempty( terr ),
        temp = { [ tsolv.name, ': ', terr ] };
    end
    reject = [ reject, temp(:)' ]; %#ok
end
if ~isempty( reject ),
    reject = sprintf( '    %s\n', reject{:} );
    warning( 'CVX:SolverErrors', 'The following errors were issued when initializing the solvers:\n%sPlease check your installation and re-run CVX_SETUP.\nThese solvers are unavailable for this session.%s', reject );
end
if nrej == length( solvers ),
    clear global cvx___
    cvx_throw( 'All solvers were disabled due to various errors.\nPlease re-run CVX_SETUP and, if necessary, contact CVX Research for support.' );
elseif ~isempty(solvers(cvx___.solvers.map.default).error),
    cvx___.solvers.map.default = ndefault;
    cvx___.solvers.selected = ndefault;
    warning( 'CVX:SolverErrors', 'The default solver has temporarily been changed to %s.', solvers(ndefault).name );
end

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
