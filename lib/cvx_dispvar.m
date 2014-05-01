function [ names, sizes ] = cvx_dispvar( v, name, isdual )
[ names, sizes ] = dispvar_( v, name, isdual );
if nargout <= 1 && ~isempty( names ),
    mlen = max( cellfun( @numel, names ) ) + 3;
    spc = ' '; spc(1,2:mlen+3) = ' ';
    names = cellfun( @(x)[x,spc(1:mlen-numel(x))], names, 'UniformOutput', false );
    names = strcat( names, sizes );
end
if nargout == 0,
    names = sprintf( '   %s\n', names{:} );
    fprintf( 1, '\n%s\n', names );
    clear names
end
    
function [ names, sizes ] = dispvar_( v, name, isdual )
names = {}; 
sizes = {};
switch class( v ),
    case 'struct',
        if numel( v ) > 1,
            sv = size( v );
            sstr = struct( 'type', '()', 'subs', {cell(1,length(sv))} );
            isvec = all(sv(2:end)==1);
            for k = 1 : prod(sv),
                [ sstr.subs{:} ] = ind2sub( sv, k );
                if isvec,
                    label = sprintf( '%s(%d)', name, k );
                else
                    label = sprintf( ',%d', sstr.subs{:} );
                    label = sprintf( '%s(%s)', name, label(2:end) );
                end
                [ name2, size2 ] = dispvar_( v(k), sstr, label, isdual );
                names = [ names, name2 ]; %#ok
                sizes = [ sizes, size2 ]; %#ok
                if k == 1, name(:) = ' '; end
            end
        else
            fn = fieldnames( v );
            if ~isempty( name ), name( end + 1 ) = '.'; end
            names = {}; sizes = {};
            for k = 1 : length( fn ),
                [ name2, size2 ] = dispvar_( v.(fn{k}), [ name, fn{k} ], isdual );
                names = [ names, name2 ]; %#ok
                sizes = [ sizes, size2 ]; %#ok
                if k == 1, name(1:end-1) = ' '; end
            end
        end
    case 'cell',
        sv = size( v );
        sstr = cell( 1, length(sv) );
        isvec = all(sv(2:end)==1);
        for k = 1 : prod(sv),
            [ sstr{:} ] = ind2sub( sv, k );
            if isvec,
                label = sprintf( '%s{%d}', name, k );
            else
                label = sprintf( ',%d', sstr{:} );
                label = sprintf( '%s{%s}', name, label(2:end) );
            end
            [ name2, size2 ] = dispvar_( v{k}, label, isdual );
            names = [ names, name2 ]; %#ok
            sizes = [ sizes, size2 ]; %#ok
            if k == 1, name( 1 : end ) = ' '; end
        end
    case 'cvx';
        names{1} = name;
        sizes{1} = type(v);
    case 'double',
        names{1} = name;
        if isdual,
            sizes{1} = 'unassigned';
        else
            sizes{1} = [ type(cvx(v)), ' constant' ];
        end
    otherwise,
        names{1} = name;
        sizes{1} = 'unknown';
end
