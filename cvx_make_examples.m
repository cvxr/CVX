function cvx_make_examples( varargin )

%
% Determine the base path
%

odir = cd;
try,
    base = dbstack( '-completenames' );
    base = base(1);
    base = base.file;
catch,
    base = dbstack;
    base = base(1);
    base = base.name;
end
base = fullfile( fileparts( base ), 'examples' );

%
% Check the force flag
%

args = varargin;
temp = strcmp( args, '-force' );
force = any( temp );
if force, args(temp) = []; end
if isempty( args ), args = { base }; end

%
% Process the arguments
%

for k = 1 : length( args ),
    
    %
    % Check the validity of the file or directory
    %
    
    file = args{k};
    switch exist( file ),
        case 0,
            error( sprintf( 'Cannot find file or directory: %s', file ) );
        case 2,
            [ mpath, file, ext, versn ] = fileparts( which( file ) );
            file = [ file, ext, versn ];
            if ~strcmp( ext, '.m' ),
                error( sprintf( 'Must be an m-file: %s', file ) );
            end
        case 7,
            cd( file );
            mpath = cd;
            cd( odir );
            file = '';
        otherwise,
            error( sprintf( 'Invalid file: %s', file ) );
    end
    if length( mpath ) < length( base ) | strncmpi( mpath, base, length( base ) ) == 0,
        error( sprintf( 'Not a valid a subdirectory of cvx/examples/: %s', mpath ) );
    end
    
    %
    % Process the file or directory
    %   

    try,
        if isempty( file ),
            generate_directory( mpath, '', force );
        else,
            cd( mpath );
            generate_file( file, '', force );
        end
        cd( odir );
    catch,
        err = lasterror;
        cd( odir );
        rethrow( err );
    end
    
end

function title = generate_directory( mpath, prefix, force )

disp( sprintf( '%sDirectory: %s', prefix, mpath ) );
prefix = [ prefix, '   ' ];
cd( mpath );
mpath = cd;
names = dir;
mlen  = 0;
isdir = [ names.isdir ];
dates = { names.date };
names = { names.name };
files = struct( 'name', {}, 'title', {} );
for k = 1 : length( names ),
    name = names{k};
    cd( mpath );
    if isdir( k ),
        if name(1) == '.', continue; end
        if strcmp( name, 'html' ), continue; end
        temp = generate_directory( name, prefix, force );
        name( end + 1 ) = '/';
    else,
        if length( name ) < 2 | ~strcmp( name(end-1:end), '.m' ), continue; end
        if strcmp( name, 'Contents.m' ), continue; end
        temp = generate_file( name, prefix, force );
    end
    if isempty( temp ), temp = '(no title)'; end
    files( end + 1 ) = struct( 'name', name, 'title', temp );
    mlen = max( mlen, length( name ) );
end
cd( mpath );
title = '';
comments = {};
if isempty( dir( 'Contents.m' ) ),
    fidr = -1;
else,
    fidr = fopen( 'Contents.m', 'r' );
end
if fidr >= 0 | ~isempty( files ) | ~isempty( dirs ),
    fprintf( 1, '%sContents.m ...', prefix );
    if fidr >= 0,
        title = fgetl( fidr );
        if length( title ) > 2 & title( 1 ) == '%' & title( 2 ) == ' ' & title( 3 ) ~= ' ',
            title = title( min( find( title ~= '%' & title ~= ' ' ) ) : end );
        else,
            title = '';
        end
    else,
        title = '';
    end
    [ fidw, message ] = fopen( 'Contents.m.new', 'w+' );
    if fidw < 0,
        if fidr >= 0, fclose( fidr ); end
        error( 'Cannot open Contents.m.new\n   %s', message );
    end
    if ~isempty( title ),
        fprintf( fidw, '%% %s\n', title );
        while ~feof( fidr ),
            temp = fgetl( fidr );
            if temp( 1 ) ~= '%' | ~any( temp ~= '%' & temp ~= ' ' ), break; end
            fprintf( fidw, '%s\n', temp );
        end
        fprintf( fidw, '%%\n' );
    end
    for k = 1 : length( files ),
        tfile = files(k);
        tfile.name(end+1:mlen) = ' ';
        if isempty( tfile.title ),
            fprintf( fidw, '%%  %s - (no title)\n', tfile.name );
        else,
            fprintf( fidw, '%%  %s - %s\n', tfile.name, tfile.title );
        end
    end
    fprintf( fidw, 'help Contents\n' );
    fseek( fidw, 0, -1 );
    c2 = fread( fidw, Inf, 'uint8' );
    fclose( fidw );
    if fidr >= 0,
        fseek( fidr, 0, -1 );
        c1 = fread( fidr, Inf, 'uint8' );
        fclose( fidr );
    else,
        c1 = [];
    end
    if length( c1 ) ~= length( c2 ) | any( c1 ~= c2 ),
        [ success, message ] = movefile( 'Contents.m.new', 'Contents.m', 'f' );
        if ~success,
            error( 'Cannot move Contents.m.new into place\n   %s', message );
            delete Contents.m.new
        end
        if fidr >= 0,
            fprintf( 1, ' updated.\n' );
        else, 
            fprintf( 1, ' created.\n' );
        end
    else,
        delete Contents.m.new
        fprintf( 1, ' up to date.\n' );
    end
end

function title = generate_file( name, prefix, force )

if length( name ) < 2 | ~strcmp( name(end-1:end), '.m' ),
    error( 'Not an m-file.' );
elseif strcmp(name(end-1:end), 'Contents.m' ),
    error( 'To generate the Contents.m file, you must run this function on the entire directory.' );
else,
    fprintf( 1, '%s%s: ', prefix, name );
end

dd = dir( name );
ndate = date_convert( dd.date );
[ fidr, message ] = fopen( name, 'r' );
if fidr < 0,
    error( 'Cannot open the source file\n   %s', message );
end
title = fgetl( fidr );
if title( 1 ) ~= '%',
    fclose( fidr );
    fprintf( 1, 'not in publishable form (perhaps a function?).\n' );
    title = '';
    return
else,
    title = title( min( find( title ~= '%' & title ~= ' ' ) ) : end );
end
hfile = [ name(1:end-1), 'html' ];
odir = cd;
hdir = 'html';
hdate = 0;
if exist( hdir ) ==  7,
    cd( hdir );
    df = dir( hfile );
    if length( df ) == 1,
        hdate = date_convert( df.date );
    end
    cd( odir );
end
if force | hdate <= ndate,
    if hdate == 0,
        fprintf( 1, 'creating %s ...', hfile );
    else,
        fprintf( 1, 'updating %s ...', hfile );
    end
    try,
        name2 = [ name, '__' ];
        [ fidw, message ] = fopen( name2, 'w+' );
        if fidw < 0,
            error( 'Cannot open the temporary file\n   %s', message );
        end
        fprintf( fidw, '%%%% %s\n\n', title );
        fwrite( fidw, fread( fidr, Inf, 'uint8' ), 'uint8' );
        fclose( fidw );
        cvx_clear
        publish( name2, struct( 'format', 'html', 'useNewFigure', 0, 'stopOnError', 1, 'createThumbnail', 0 ) );
        delete( name2 );
        if ~isempty( dir( name2 ) ),
            error( 'Cannot delete the temporary file' );
        end
        fclose( fidr );
        cd( odir );
        fprintf( 1, ' done.\n' );
    catch,
        err = lasterror;
        fprintf( 1, ' aborted.\n' );
        fclose( fidr );
        cd( odir );
        rethrow( err );
    end
    close all
else,
    fprintf( 1, 'up to date.\n' );
end

function dnum = date_convert( dstr )
persistent mstrs
if isempty( mstrs ),
    mstrs = { 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' };
end
% DD-MMM-YY HH:MM:SS
S = sscanf( dstr, '%d-%3s-%d %d:%d:%d' );
S = [ S(5), find(strcmp(char(S(2:4)'),mstrs)), S(1), S(6), S(7), S(8) ];
dnum = S(6) + 100 * ( S(5) + 100 * ( S(4) + 100 * ( S(3) + 100 * ( S(2) + 100 * S(1) ) ) ) );
