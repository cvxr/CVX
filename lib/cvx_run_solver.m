function varargout = cvx_run_solver( sfunc, varargin )
settings_arg = varargin{end};
settings     = varargin{end-1};
inputs       = varargin(1:end-nargout-2);
dumpfile = '';
custom_on = false;
if isstruct( settings ),
    for f = fieldnames( settings )',
        sval = settings.(f{1});
        if isequal( f{1}, 'dumpfile' ),
            dumpfile = sval;
        else
            custom_on = true;
            inputs{settings_arg}.(f{1}) = sval;
        end
    end
end
if custom_on,
    fprintf( 'NOTE: custom settings have been set for this solver.\n' );
end
if ~isempty( dumpfile ),
    if ~ischar( dumpfile ) || size( dumpfile, 1 ) > 1,
        error( 'CVX:Dumpfile', 'Invalid filename for the dumpfile.' );
    elseif length(dumpfile) < 4 || ~strcmpi(dumpfile(end-3:end),'.mat'),
        dumpfile = [ dumpfile, '.mat' ];
    end
    fprintf( 'Saving output to: %s\n', dumpfile );
    fprintf( '------------------------------------------------------------\n');
    otp_names = varargin(end-nargout-1:end-2);
    inp_names = cell(1,length(inputs));
    for k = 1 : length(inp_names),
        inp_names{1,k} = inputname(k+1);
    end
    fid = fopen( dumpfile, 'w+' );
    if fid == 0,
        error( 'CVX:Dumpfile', 'Cannot open file %s for writing\n', fid );
    end
    fclose( fid );
    diary( dumpfile );
elseif custom_on,
    fprintf( '------------------------------------------------------------\n');
end
errmsg = [];
try
   [ varargout{1:nargout} ] = sfunc( inputs{:} );
catch errmsg
   [ varargout{1:nargout} ] = deal( [] );
end
if ~isempty( dumpfile ),
    diary( 'off' );
    fid = fopen( dumpfile, 'r' );
    if fid ~= 0,
        output = fread( fid, Inf, '*char' )';
    else
        output = '<Could not save>';
    end
    dstruct = cell2struct( [ inputs, varargout, output ], [ inp_names, otp_names, 'output' ], 2 ); %#ok
    save( dumpfile, '-struct', 'dstruct' );
end
if ~isempty( errmsg ),
    rethrow( errmsg );
end
