function args = cvx_parse
try
    global cvx___ %#ok
    oargs = cvx___.args{1};
    nargs = length( oargs );
    if isequal( cvx___.args{2}, 'dual' ),
        nvars = nargs;
        rpat = '^([a-zA-Z]\w*)({.*})?$';
        darg = {};
    else
        nvars = cvx___.args{2};
        darg = cvx___.args{3};
        rpat = '^([a-zA-Z]\w*)(\(.*\))?$';
        if nvars < 0, nvars = nargs; end
    end
    args = struct( 'name', cell( 1, nargs ), 'args', cell( 1, nargs ) );
    for k = 1 : nargs,
        tok = regexp( oargs{k}, rpat, 'tokens' );
        if isempty( tok ),
            if k <= nvars, type = 'variable'; else type = 'structure'; end
            cvx_throw( 'Invalid %s specification: %s', type, oargs{k} );
        end
        tok = tok{1};
        nam = tok{1};
        if length(cvx___.problems) > 1 && ( nvars < 0 || k <= nvars ),
            if ~isvarname( nam ),
                cvx_throw( 'Invalid variable name: %s', nam );
            elseif nam(end) == '_',
                cvx_throw( 'Invalid variable name: %s\n   Variables ending in underscores are reserved for internal use.', nam );
            elseif isfield( cvx___.reswords, nam ),
                if cvx___.reswords.(nam) == 'S',
                    cvx_throw( 'Invalid variable name: %s\n   This is a reserved word in CVX.\n   Trying to declare a structured matrix? Use the VARIABLE keyword instead.', nam );
                else
                    cvx_throw( 'Invalid variable name: %s\n   This is a reserved word in CVX.', nam );
                end
            end
            switch exist( nam ), %#ok
                case {0,1},
                case 5,
                    cvx_throw( 'Variable name "%s" is the name of a built-in MATLAB function.\nPlease choose a different name.', nam );
                case 8,
                    cvx_throw( 'Variable name "%s" is the name of an existing MATLAB class.\nPlease choose a different name.', nam );
                otherwise,
                    mpath = which( nam );
                    if ~isempty( nam ),
                        if strncmp( nam, matlabroot, length(matlabroot) ),
                            cvx_throw( 'Variable name "%s" is the name of an existing MATLAB function or directory:\n    %s\nPlease choose a different name.', nam, mpath );
                        elseif strncmp( mpath, cvx___.where, length(cvx___.where) ),
                            cvx_throw( 'Variable name "%s" matches the name of a CVX function or directory:\n    %s\nPlease choose a different name.', nam, mpath );
                        else
                            warning( 'Variable name "%s" matches the name of an function or directory:\n    %s\nThis may cause unintended behavior with CVX models.\nPlease consider moving this file or choosing a different variable name.', nam, mpath );
                        end
                    end
            end
        end
        if length(tok) < 2 || isempty( tok{2} ),
            if k <= nvars, arg = darg;
            else arg = {}; end
        else
            try
                if k <= nvars,
                    arg = evalin( 'caller', [ '[', tok{2}(2:end-1), ']' ] );
                    arg = cvx_get_dimlist( arg );
                else
                    arg = evalin( 'caller', [ '{', tok{2}(2:end-1), '}' ] );
                end
            catch exc
                cvx_throw( 'Error attempting to evaluate arguments of: %s\n   %s', oargs{k}, exc.message );
            end
        end
        args(k).name = nam;
        args(k).args = arg;
        args(k).orig = oargs{k};
    end
    cvx___.args = {};
catch exc
    cvx___.args = {};
    rethrow( exc );
end
