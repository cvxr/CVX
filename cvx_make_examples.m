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
    file = args{k};
    if any( file == '*' ),
        files = dir( file );
        files = { files.name };
    else,
        files = { file };
    end
    for j = 1 : length( files );
        
        %
        % Check the validity of the file or directory
        %

        file = files{j};
        switch exist( file ),
            case 0,
                error( sprintf( 'Cannot find file or directory: %s', file ) );
            case 2,
                [ mpath, file, ext, versn ] = fileparts( which( file ) );
                file = [ file, ext, versn ];
                if ~strcmp( ext, '.m' ),
                    error( sprintf( 'Must be an m-file: %s', file ) );
                elseif strcmp( file, 'Contents.m' ) & length( files ) > 1,
                    continue;
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

        if isempty( file ) & strcmp( mpath, base ),
            cd( base );
            [ fidr, message ] = fopen( 'index.html', 'r' );
            if fidr < 0,
                error( 'Cannot open index.html\n   %s', message );
            end
            t_head = fread( fidr, Inf, 'uint8=>char' )';
            t_tail = strfind( t_head, '<!-- END TOC -->' );
            if length( t_tail ) ~= 1
                error( 'Corrupt index.html file, cannot proceed.' );
            end
            t_tail = t_head( t_tail : end );
            n_head = strfind( t_head, '<!-- BEGIN TOC -->' );
            if length( n_head ) ~= 1,
                error( 'Corrupt index.html file, cannot proceed.' );
            end
            t_head( n_head + length( '<!-- BEGIN TOC -->' ) + 1 : end ) = [];
            fclose( fidr );
            [ fidw, message ] = fopen( 'index.html.new', 'w+' );
            if fidw < 0,
                error( 'Cannot open index.html.new\n   %s', message );
            end
            fwrite( fidw, t_head, 'char' );
            fprintf( fidw, '<h3>Last updated: %s</h3>\n', date );
            fprintf( fidw, '<a href="#" onClick="expandTree(''tree1''); return false;">Expand all</a>&nbsp;&nbsp;\n' );
            fprintf( fidw, '<a href="#" onClick="collapseTree(''tree1''); return false;">Collapse all</a>\n' );
        else,
            fidw = -1;
        end

        if isempty( file ),
            generate_directory( mpath, '', force, fidw, base );
        else,
            cd( mpath );
            generate_file( file, '', force );
        end
        cd( odir );

        if fidw >= 0,
            fwrite( fidw, t_tail, 'char' );
            fclose( fidw );
            cd( mpath )
            compare_and_replace( '', 'index.html' );
        end
    end
end

function [ title, files ] = generate_directory( mpath, prefix, force, fidc, base )

disp( sprintf( '%sDirectory: %s', prefix, mpath ) );
prefix = [ prefix, '   ' ];
cd( mpath );
mpath = cd;

%
% Open Contents.m file and retrieve title and comments
%

title = '';
comments = {};
[ fidr, message ] = fopen( 'Contents.m', 'r' );
if fidr >= 0,
    temp = fgetl( fidr );
    if length( temp ) > 2 & temp( 1 ) == '%' & temp( 2 ) == ' ' & temp( 3 ) ~= ' ',
        title = temp( min( find( temp ~= '%' & temp ~= ' ' ) ) : end );
        while ~feof( fidr ),
            temp = fgetl( fidr );
            if temp( 1 ) ~= '%' | ~any( temp ~= '%' & temp ~= ' ' ), break; end
            comments{end+1} = temp( min( find( temp ~= '%' & temp ~= ' ' ) ) : end );
        end
    end
    fclose( fidr );
elseif ~isempty( dir( 'Contents.m' ) ),
    error( 'Cannot open Contents.m for reading\n   %s', message );
end

%
% Open Comments.m.new and deposit title and comments
%

[ fidw, message ] = fopen( 'Contents.m.new', 'w+' );
if fidw < 0,
    if fidr >= 0, fclose( fidr ); end
    error( 'Cannot open Contents.m.new\n   %s', message );
elseif ~isempty( title ),
    fprintf( fidw, '%% %s\n', title );
    for k = 1 : length( comments ),
        fprintf( fidw, '%% %s\n', comments{k} );
    end
    fprintf( fidw, '%%\n' );
end

%
% Determine the names of the files and subdirectories to process
%

dd = dir;
dnames = {};
fnames = {};
files = struct( 'name', {}, 'title', {} );
mlen = 0;
for k = 1 : length( dd ),
    name = dd(k).name;
    if dd(k).isdir,
        if name(1) == '.' | strcmp( name, 'html' ), continue; end
        dnames{end+1} = name;
    else,
        if length( name ) < 2 | ~strcmp( name(end-1:end), '.m' ), continue; end
        if strcmp( name, 'Contents.m' ), continue; end
        fnames{end+1} = name;
    end
end

%
% If needed, deposit directory name, title, and comments into index.html
%

if fidc >= 0,
    dpath = mpath( length(base) + 2 : end );
    dpath(dpath=='\') = '/';
    if ~isempty( dpath ),
        dpath(end+1) = '/';
        [ dpath2, dname ] = fileparts( mpath );
        if isempty( dnames ),
            mclass = 'liClosed';
        else,
            mclass = 'liOpen';
        end
        if isempty( title ),
            fprintf( fidc, '<li class="%s"><a href="%s">%s/</a>\n', mclass, dpath, dname ),
        else,
            fprintf( fidc, '<li class="%s"><b>%s</b>\n', mclass, title );
%            fprintf( fidc, '<li class="%s"><b>%s</b> (<a href="%s">%s/</a>)\n', mclass, title, dpath, dname );
        end
        if false,
            if ~isempty( comments ),
                fprintf( fidc, '<br />' );
                for k = 1 : length( comments ),
                    fprintf( fidc, '%s\n', comments{k} );
                end
            end
        end
        fprintf( fidc, '<ul>\n' );
    else,
        fprintf( fidc, '<ul class="mktree" id="tree1">\n' );
    end
end

%
% Process the subdirectories
%

for k = 1 : length( dnames ),
    name = dnames{k};
    if name(1) == '.' | strcmp( name, 'html' ), continue; end
    [ temp, contents ] = generate_directory( name, prefix, force, fidc, base );
    if isempty( temp ), temp = '(no title)'; end
    name( end + 1 ) = '/';
    files( end + 1 ) = struct( 'name', name, 'title', temp );
    mlen = max( mlen, length( name ) );
    cd( mpath );
end

%
% Process the files
%

if fidc >= 0 & isempty( dpath ) & ~isempty( fnames ) & ~isempty( dnames ),
    need_misc = true;
    fprintf( fidc, '<li><b>Uncategorized files</b><ul>\n' );
else,
    need_misc = false;
end
for k = 1 : length( fnames ),
    name = fnames{k};
    if length( name ) < 2 | ~strcmp( name(end-1:end), '.m' ), continue; end
    if strcmp( name, 'Contents.m' ), continue; end
    temp = generate_file( name, prefix, force );
    if fidc >= 0,
        if isempty( temp ),
            fprintf( fidc, '<li><a href="%s%s">%s</a></li>\n', dpath, name, name );
        else,
            fprintf( fidc, '<li><a href="%shtml/%shtml">%s</a> (<a href="%s%s">%s</a>)</li>\n', dpath, name(1:end-1), temp, dpath, name, name );
        end
    end
    if isempty( temp ), temp = '(no title)'; end
    files( end + 1 ) = struct( 'name', name, 'title', temp );
    mlen = max( mlen, length( name ) );
end
if need_misc,
    fprintf( fidc, '</ul></li>\n' );
elseif isempty( fnames ) & isempty( dnames ),
    fprintf( fidc, '<li>(no files)</li>\n' );
end

%
% Finish catalog.html.new entry
%

if fidc >= 0,
    fprintf( fidc, '</ul>\n' );
    if ~isempty( dpath ),
        fprintf( fidc, '</li>\n' );
    end
end

%
% Finish Contents.m.new
%

for k = 1 : length( files ),
    tfile = files(k);
    tfile.name(end+1:mlen) = ' ';
    if isempty( tfile.title ),
        fprintf( fidw, '%%  %s - (no title)\n', tfile.name );
    else,
        fprintf( fidw, '%%  %s - %s\n', tfile.name, tfile.title );
    end
end
fclose( fidw );

%
% Compare Contents.m and Contents.m.new and update if necessary
%

cd( mpath )
compare_and_replace( prefix, 'Contents.m' );

function title = generate_file( name, prefix, force )

if length( name ) < 2 | ~strcmp( name(end-1:end), '.m' ),
    error( 'Not an m-file.' );
elseif strcmp( name, 'Contents.m' ),
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
        if isempty( title ),
            fprintf( fidw, '%%%% %s\n\n', name );
        else,
            fprintf( fidw, '%%%% %s\n\n', title );
            while ~feof( fidr ),
                temp = fgetl( fidr );
                if temp( 1 ) ~= '%' | any( temp( 2 : end ) ~= ' ' ),
                    fprintf( fidw, '%s\n', temp );
                    break;
                end
            end
        end
        fwrite( fidw, fread( fidr, Inf, 'uint8' ), 'uint8' );
        fclose( fidw );
        cvx_clear
        cvx_quiet( false );
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

function compare_and_replace( prefix, oldname )

names = { oldname, [ oldname, '.new' ] };
fprintf( 1, '%s%s ... ', prefix, oldname );
fids = [];
c = {};
for k = 1 : 2,
    [ fids(k), message ] = fopen( names{k}, 'r' );
    if fids(k) < 0 & ~isempty( dir( names{k} ) ),
        error( 'Cannot open file %s for reading:\n   %s', names{k}, message );
    end
    c{k} = fread( fids(k), Inf, 'uint8' );
    fclose( fids(k) );
end
if length( c{2} ) == 0,
    if fids(k) >= 0,
        fprintf( 1, ' removed.\n' );
        delete( oldname );
    end
    delete( names{2} );
elseif length( c{1} ) ~= length( c{2} ) | any( c{1} ~= c{2} ),
    [ success, message ] = movefile( names{2}, names{1}, 'f' );
    if ~success,
        error( 'Cannot move %s into place\n   %s', names{2}, message );
        delete( names{2} )
    end
    if ~isempty( c{1} ),
        fprintf( 1, ' updated.\n' );
    else,
        fprintf( 1, ' created.\n' );
    end
else,
    delete( names{2} )
    fprintf( 1, ' up to date.\n' );
end
