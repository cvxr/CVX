function [ prob, name, args ] = cvx_parse
try
    global cvx___ %#ok
    prob = evalin( 'caller', 'cvx_problem', '[]' );
    if ~isa( prob, 'cvxprob' ),
        error( 'CVX:NoProblem', 'No CVX model exists in this scope.' );
    end
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
    name  = cell( 1, nargs );
    args  = cell( 1, nargs );
    for k = 1 : nargs,
        tok = regexp( oargs{k}, rpat, 'tokens' );
        if isempty( tok ),
            if k <= nvars, type = 'Variable'; else type = 'Structure'; end
            error( sprintf('CVX:%s',type), 'Invalid %s specification: %s', lower(type), oargs{k} );
        end
        tok = tok{1};
        nam = tok{1};
        if nvars < 0 || k <= nvars,
            if ~isvarname( nam ),
                error( 'CVX:Variable', 'Invalid variable name: %s', nam );
            elseif nam(end) == '_',
                error( 'CVX:Variable', 'Invalid variable name: %s\n   Variables ending in underscores are reserved for internal use.', nam );
            elseif isfield( cvx___.reswords, nam ),
                if cvx___.reswords.(nam) == 'S',
                    error( 'CVX:Variable', 'Invalid variable name: %s\n   This is a reserved word in CVX.\n   Trying to declare a structured matrix? Use the VARIABLE keyword instead.', nam );
                else
                    error( 'CVX:Variable', 'Invalid variable name: %s\n   This is a reserved word in CVX.', nam );
                end
            end
            switch evalin('caller',['exist(''',nam,''')']),
                case {0,1},
                case 5,
                    error( 'CVX:Variable', 'Variable name "%s" is the name of a built-in MATLAB function.\nPlease choose a different name.', nam );
                case 8,
                    error( 'CVX:Variable', 'Variable name "%s" is the name of an existing MATLAB class.\nPlease choose a different name.', nam );
                otherwise,
                    mpath = which( nam );
                    if ~isempty( nam ),
                        if strncmp( nam, matlabroot, length(matlabroot) ),
                            error( 'CVX:Variable', 'Variable name "%s" is the name of an existing MATLAB function or directory:\n    %s\nPlease choose a different name.', nam, mpath );
                        elseif strncmp( mpath, cvx___.where, length(cvx___.where) ),
                            error( 'CVX:Variable', 'Variable name "%s" matches the name of a CVX function or directory:\n    %s\nPlease choose a different name.', nam, mpath );
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
                    arg = cvx_get_dimlist( arg, 'default', [1,1] );
                else
                    arg = evalin( 'caller', [ '{', tok{2}(2:end-1), '}' ] );
                end
            catch exc
                error( exc.identifier, 'Error attempting to evaluate arguments of: %s\n   %s', oargs{k}, exc.message );
            end
        end
        name{k} = nam;
        args{k} = arg;
    end
    cvx___.args = {};
catch exc
    cvx___.args = {};
    rethrow( exc );
end
