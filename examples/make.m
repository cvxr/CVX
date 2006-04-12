function make( varargin )

%
% Quick wrapper to save base workspace
%

evalin( 'base', 'global cvx___ws___' );
vars = evalin( 'base', 'who' );
varg = evalin( 'base', 'who( ''global'' )' );
varg = intersect( varg, vars );
vars = setdiff( vars, varg );
temp = sprintf( '%s,', vars{:} );
evalin( 'base', [ 'cvx___ws___ = {', temp(1:end-1), '};' ] );
evalin( 'base', 'clear' );
try
    make2( varargin{:} );
    lerr = [];
catch
    lerr = lasterror;
end
evalin( 'base', 'clear' );
if ~isempty( varg ),
    evalin( 'base', sprintf( '%s ', 'global', varg{:} ) );
end
if ~isempty( vars ),
    evalin( 'base', [ '[', temp(1:end-1), '] = deal( cvx___ws___{:} );' ] );
end
clear global cvx___ws___
if ~isempty( lerr ),
    rethrow( lerr );
end

function make2( varargin )

%
% Determine the base path
%

odir = cd;
try
    base = dbstack( '-completenames' );
    base = base(1);
    base = base.file;
catch
    base = dbstack;
    base = base(1);
    base = base.name;
end
base = fileparts( base );
fclose all;
close all;

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
    else
        files = { file };
    end
    for j = 1 : length( files );

        %
        % Check the validity of the file or directory
        %

        file = files{j};
        switch exist( file, 'file' ),
            case 0,
                error( 'Cannot find file or directory: %s', file );
            case 2,
                [ mpath, file, ext, versn ] = fileparts( which( file ) );
                file = [ file, ext, versn ];
                if ~strcmp( ext, '.m' ),
                    error( 'Must be an m-file: %s', file );
                elseif strcmp( file, 'Contents.m' ) && length( files ) > 1,
                    continue;
                elseif strcmp( file, 'make.m' ) && strcmp( mpath, base ),
                    continue;
                end
            case 7,
                cd( file );
                mpath = cd;
                cd( odir );
                file = '';
            otherwise,
                error( 'Invalid file: %s', file );
        end
        if length( mpath ) < length( base ) || strncmpi( mpath, base, length( base ) ) == 0,
            error( 'Not a valid a subdirectory of cvx/examples/: %s', mpath );
        end

        %
        % Process the file or directory
        %

        if isempty( file ) && strcmp( mpath, base ),
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
            fprintf( fidw, '<p><a href="#" onclick="expandTree(''tree1''); return false;">Expand all</a>&nbsp;&nbsp;\n' );
            fprintf( fidw, '<a href="#" onclick="collapseTree(''tree1''); return false;">Collapse all</a></p>\n' );
        else
            fidw = -1;
        end

        if isempty( file ),
            generate_directory( mpath, '', force, fidw, base );
        else
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

persistent htmlsrc htmldst
if isempty( htmlsrc ),
    htmlsrc = { '&', '<', '>', 'http://www.stanford.edu/([^\s\),]*)' };
    htmldst = { '&amp;', '&lt;', '&gt;', '<a href="http://www.stanford.edu/$1">$1</a>' };
end

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
    if length( temp ) > 2 && temp( 1 ) == '%' && temp( 2 ) == ' ' && temp( 3 ) ~= ' ',
        title = temp( min( find( temp ~= '%' & temp ~= ' ' ) ) : end );
        while ~feof( fidr ),
            temp = fgetl( fidr );
            if temp( 1 ) ~= '%' || ~any( temp ~= '%' & temp ~= ' ' ), break; end
            comments{end+1} = temp( min( find( temp ~= '%' & temp ~= ' ' ) ) : end );
        end
    end
    fclose( fidr );
elseif ~isempty( dir( 'Contents.m' ) ),
    error( 'Cannot open Contents.m for reading\n   %s', message );
end

%
% Read the entries, and process the scripts and functions
%

dd = dir;
mlen = 0;
files = struct( 'name', {}, 'title', {}, 'type', {} );
for k = 1 : length( dd ),
    name = dd(k).name;
    if dd(k).isdir,
        if name(1) == '.' || strcmp( name, 'html' ), continue; end
        name(end+1) = '/';
        files( end + 1 ) = struct( 'name', name, 'title', '', 'type', 'dir' );
    elseif length( name ) > 2,
        ndx = max(find(name=='.'));
        if isempty( ndx ), continue; end
        switch name(ndx+1:end),
            case 'm',
                if strcmp( name, 'Contents.m' ) || strcmp( name, 'make.m' ), continue; end
                [ temp, isfunc ] = generate_file( name, prefix, force );
                if isfunc, type = 'func'; else type = 'script'; end
                files( end + 1 ) = struct( 'name', name, 'title', temp, 'type', type );
            case 'tex',
                temp = generate_doc( name, prefix, force );
                files( end + 1 ) = struct( 'name', name, 'title', temp, 'type', 'tex' );
            case { 'pdf', 'ps' },
                if any( strcmp( { dd.name }, [name(1:ndx+1),'tex'] ) ), continue; end
                files( end + 1 ) = struct( 'name', name, 'title', '', 'type', 'doc' );
            otherwise,
                continue;
        end
    end
    mlen = max( mlen, length(name) );
end

%
% Sort the files
%

if ~isempty( files ),
    [ fnames, ndxs ] = sort( { files.title } );
    files = files(ndxs);
    ftypes = { files.type };
    tdir  = strcmp( ftypes, 'dir' );
    tfun  = strcmp( ftypes, 'func' );
    tdoc  = strcmp( ftypes, 'doc' ) | strcmp( ftypes, 'tex' );
    tscr  = ~( tdir | tfun | tdoc );
    t1    = strncmp( fnames, 'Exercise', 8 ) & tscr;
    t2    = strncmp( fnames, 'Example',  7 ) & tscr;
    t3    = strncmp( fnames, 'Section',  7 ) & tscr;
    t4    = strncmp( fnames, 'Figure', 6 ) & tscr;
    t5    = ~( t1 | t2 | t3 | t4 ) & tscr;
    ndxs  = [ find(tdir(:)); find(t3(:)); find(t2(:)); find(t4(:)); find(t5(:)); find(t1(:)) ; find(tfun(:)) ; find(tdoc(:)) ];
    files = files(ndxs);
end

%
% Fill out the index.html file
%

if fidc >= 0,

    dpath = mpath( length(base) + 2 : end );
    dpath(dpath=='\') = '/';
    if ~isempty( dpath ),
        dpath(end+1) = '/';
        [ dpath2, dname ] = fileparts( mpath );
        if ~isempty( files ) && any( tdir ),
            mclass = 'liOpen';
        else
            mclass = 'liClosed';
        end
        if isempty( title ),
            fprintf( fidc, '<li class="%s"><a href="%s">%s/</a>\n', mclass, dpath, dname ),
        else
            fprintf( fidc, '<li class="%s"><b>%s</b>\n', mclass, regexprep( title, htmlsrc, htmldst ) );
        end
        if ~isempty( comments ),
            fprintf( fidc, '<blockquote>\n' );
            for k = 1 : length( comments ),
                fprintf( fidc, '%s\n', regexprep( comments{k}, htmlsrc, htmldst ) );
            end
            fprintf( fidc, '</blockquote>\n' );
        end
        fprintf( fidc, '<ul>\n' );
    else
        fprintf( fidc, '<ul class="mktree" id="tree1">\n' );
    end

    if isempty( files ),

        fprintf( fidc, '<li>(no files)</li>\n' );

    else

        if any( tdir ),
            for k = 1 : length( files ),
                if strcmp( files(k).type, 'dir' ),
                    files(k).title = generate_directory( files(k).name(1:end-1), prefix, force, fidc, base );
                    cd(mpath);
                end
            end
        end

        if any( tscr ),
            need_misc = isempty( dpath ) & any( tdir );
            if need_misc,
                fprintf( fidc, '<li><b>Miscellaneous examples</b><ul>\n' );
            end
            for k = 1 : length( files ),
                if strcmp( files(k).type, 'script' ),
                    name = files( k ).name;
                    temp = files( k ).title;
                    if isempty( temp ),
                        fprintf( fidc, '<li><a href="%s%s">%s</a></li>\n', dpath, name, name );
                    else
                        fprintf( fidc, '<li><a href="%shtml/%shtml">%s</a> (<a href="%s%s">%s</a>)</li>\n', dpath, name(1:end-1), regexprep( temp, htmlsrc, htmldst ), dpath, name, name );
                    end
                end
            end
            if need_misc,
                fprintf( fidc, '</ul></li>\n' );
            end
        end

        if any( tfun ),
            fprintf( fidc, '<li class="liOpen">Utility functions:<ul>\n' );
            for k = 1 : length( files ),
                if strcmp( files(k).type, 'func' ),
                    name = files( k ).name;
                    temp = files( k ).title;
                    if isempty( temp ),
                        fprintf( fidc, '<li><a href="%s%s">%s</a></li>\n', dpath, name, name );
                    else
                        fprintf( fidc, '<li><a href="%s%s">%s (%s)</a></li>\n', dpath, name, regexprep( temp, htmlsrc, htmldst ), name );
                    end
                end
            end
            fprintf( fidc, '</ul></li>\n' );
        end

        if any( tdoc ),
            fprintf( fidc, '<li class="liOpen">Documentation:<ul>\n' );
            for k = 1 : length( files ),
                if strcmp( files(k).type, 'doc' ) || strcmp( files(k).type, 'tex' ),
                    name = files( k ).name;
                    if strcmp( files(k).type, 'tex' ),
                        name = [ name(1:end-4), 'pdf' ];
                    end
                    temp = files( k ).title;
                    if isempty( temp ),
                        fprintf( fidc, '<li><a href="%s%s">%s</a></li>\n', dpath, name, name );
                    else
                        fprintf( fidc, '<li><a href="%s%s">%s (%s)</a></li>\n', dpath, name, regexprep( temp, htmlsrc, htmldst ), name );
                    end
                end
            end
            fprintf( fidc, '</ul></li>\n' );
        end

    end

    fprintf( fidc, '</ul>\n' );
    if ~isempty( dpath ),
        fprintf( fidc, '</li>\n' );
    end

end

%
% Create Contents.m.new
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
for k = 1 : length( files ),
    tfile = files(k);
    tfile.name(end+1:mlen) = ' ';
    if isempty( tfile.title ),
        fprintf( fidw, '%%  %s - (no title)\n', tfile.name );
    else
        fprintf( fidw, '%%  %s - %s\n', tfile.name, tfile.title );
    end
end
fprintf( fidw, 'help Contents\n' );
fclose( fidw );

%
% Compare Contents.m and Contents.m.new and update if necessary
%

cd( mpath )
compare_and_replace( prefix, 'Contents.m' );

function [ title, isfunc ] = generate_file( name, prefix, force )

if length( name ) < 2 || ~strcmp( name(end-1:end), '.m' ),
    error( 'Not an m-file.' );
elseif strcmp( name, 'Contents.m' ),
    error( 'To generate the Contents.m file, you must run this function on the entire directory.' );
else
    fprintf( 1, '%s%s: ', prefix, name );
end

dd = dir( name );
ndate = date_convert( dd.date );
[ fidr, message ] = fopen( name, 'r' );
if fidr < 0,
    error( 'Cannot open the source file\n   %s', message );
end
title = '';
isfunc = false;
lasttitle = false;
founddata = false;
prefixes = {};
while ~feof( fidr ) && ( ~founddata || isempty( title ) || lasttitle ),
    temp1 = fgetl( fidr );
    temp2 = strtrim( temp1 );
    if isempty( temp2 ),
        if lasttitle, continue; end
    elseif temp2( 1 ) == '%',
        temp3 = strtrim( temp2( min( find( temp2 ~= '%' ) ) : end ) );
        if isempty( temp3 ),
            if lasttitle, continue; end
        elseif isempty( title ),
            title = temp3;
            lasttitle = true;
            continue;
        else
            lasttitle = false;
        end
    else
        lasttitle = false;
        founddata = true;
        if strncmp( temp2, 'function', 8 ) && ( length( temp2 ) == 8 || ~isvarname( temp2( 1 : 9 ) ) ),
            isfunc = true;
        end
    end
    prefixes{end+1} = temp1;
end
if isfunc,
    fprintf( 1, 'function.\n' );
    fclose( fidr );
    return
end
hfile = [ name(1:end-1), 'html' ];
odir = cd;
hdir = 'html';
hdate = 0;
if exist( hdir, 'dir' ),
    cd( hdir );
    df = dir( hfile );
    if length( df ) == 1,
        hdate = date_convert( df.date );
    end
    cd( odir );
end
if force || hdate <= ndate,
    if hdate == 0,
        fprintf( 1, 'creating %s ...', hfile );
    else
        fprintf( 1, 'updating %s ...', hfile );
    end
    name2 = [ name, '__' ];
    [ fidw, message ] = fopen( name2, 'w+' );
    if fidw < 0,
        error( 'Cannot open the temporary file\n   %s', message );
    end
    if isempty( title ),
        fprintf( fidw, '%%%% %s\n\n', name );
    else
        fprintf( fidw, '%%%% %s\n\n', title );
    end
    fprintf( fidw, '%s\n', prefixes{:} );
    fwrite( fidw, fread( fidr, Inf, 'uint8' ), 'uint8' );
    fclose( fidw );
    evalin( 'base', 'clear' );
    cvx_quiet( false );
    cvx_precision default;
    try
        publish( name2, struct( 'format', 'html', 'useNewFigure', 0, 'createThumbnail', 0 ) );
    catch
        err = lasterror;
        fprintf( 1, ' aborted.\n' );
        fclose( fidr );
        cd( odir );
        rethrow( err );
    end
    delete( name2 );
    if ~isempty( dir( name2 ) ),
        error( 'Cannot delete the temporary file' );
    end
    fclose( fidr );
    cd( odir );
    fprintf( 1, ' done.\n' );
    close all
else
    fprintf( 1, 'up to date.\n' );
end

function title = generate_doc( name, prefix, force )

if length( name ) < 5 || ~strcmp( name(end-3:end), '.tex' ),
    error( 'Not an valid TeX file.' );
else
    fprintf( 1, '%s%s: ', prefix, name );
end

dd = dir( name );
ndate = date_convert( dd.date );
[ fidr, message ] = fopen( name, 'r' );
if fidr < 0,
    error( 'Cannot open the source file\n   %s', message );
end
title = '';
while ~feof( fidr ),
    temp = strtrim( fgetl( fidr ) );
    kndx = strfind( temp, '\title{' );
    if isempty( kndx ), continue; end
    knd2 = strfind( temp(kndx(1):end), '}' );
    if isempty( knd2 ), continue; end
    title = strtrim(temp(kndx(1)+7:kndx(1)+kndx(2)-2));
    break;
end
pdffile = [ name(1:end-3), 'pdf' ];
hdate = 0;
df = dir( pdffile );
if length( df ) == 1,
    hdate = date_convert( df.date );
end
if force || hdate < ndate,
    if hdate == 0,
        fprintf( 1, 'creating %s:', hfile );
    else
        fprintf( 1, 'updating %s:', hfile );
    end
    name2 = name(1:end-4);
    eval( sprintf( '!latex %s', name2 ) );
    eval( sprintf( '!latex %s', name2 ) );
    eval( sprintf( '!bibtex %s', name2 ) );
    eval( sprintf( '!latex %s', name2 ) );
    eval( sprintf( '!latex %s', name2 ) );
    eval( sprintf( '!latex %s', name2 ) );
    eval( sprintf( '!dvips %s', name2 ) );
    eval( sprintf( '!ps2pdf %s.ps', name2 ) );
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
    if fids(k) < 0 && ~isempty( dir( names{k} ) ),
        error( 'Cannot open file %s for reading:\n   %s', names{k}, message );
    end
    c{k} = fread( fids(k), Inf, 'uint8' );
    fclose( fids(k) );
end
if isempty( c{2} ),
    if fids(k) >= 0,
        fprintf( 1, ' removed.\n' );
        delete( oldname );
    end
    delete( names{2} );
elseif length( c{1} ) ~= length( c{2} ) || any( c{1} ~= c{2} ),
    [ success, message ] = movefile( names{2}, names{1}, 'f' );
    if ~success,
        error( 'Cannot move %s into place\n   %s', names{2}, message );
        delete( names{2} )
    end
    if ~isempty( c{1} ),
        fprintf( 1, ' updated.\n' );
    else
        fprintf( 1, ' created.\n' );
    end
else
    delete( names{2} )
    fprintf( 1, ' up to date.\n' );
end
