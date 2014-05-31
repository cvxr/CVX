function [ S, itype, stypes ] = cvx_create_structure( varargin )

%CVX_CREATE_STRUCTURE Construct a basis for a structured matrix.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse the arguments, if needed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S = [];
itype = '';
stypes = {};
if isstruct( varargin{1} ),
    args = varargin{1};
    orig = { args.orig };
    name = { args.name };
    args = { args.args };
    nargs = numel( orig );
else
    nargs = nargin;
    orig = varargin;
    name = orig;
    args = cell( 1, nargs );
    if ischar( varargin{1} ),
        amin = 1;
    else
        args{1} = varargin{1};
        name{1} = '';
        varargin{1} = '';
        amin = 2;
    end
    toks = regexp( varargin, '^([a-zA-Z]\w*)(\(.*\))?$', 'tokens' );
    for k = amin : nargs,
        tok = toks{k};
        if isempty( tok ),
            if k == 1, type = 'variable'; else type = 'structure'; end
            cvx_throw( 'Invalid %s specification: %s', type, varargin{k} );
        end
        tok = tok{1};
        name{k} = tok{1};
        if length(tok) > 1 && ~isempty( tok{2} ),
            try
                args{k} = evalin( 'caller', [ '{', tok{2}(2:end-1), '};' ] );
            catch exc
                cvx_throw( 'Error evaluating structure type: %s\n    %s', orig{k}, exc.message );
            end
        else
            args{k} = {};
        end
    end
end
if nargs == 1 && nargout > 0,
    return
end
sz = args{1};
if iscell( sz ),
    sz = [ sz{:} ];
end
sz(end+1:2) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan the structure strings for symmetry and complex cases. We now handle %
% these cases specially for improved performance.                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uplo = '';
do_semi = 0;
do_skew = false;
do_conj = false;
do_comp = false;
do_symm = false;
do_nneg = false;
do_dble = false;
is_toep = false;
do_strx = false;
hflags  = false(1,nargs); n_hflags = 0;
iflags  = false(1,nargs); n_iflags = 0;
sflags  = false(1,nargs); n_sflags = 0;
bflags  = false(1,nargs); n_bflags = 0;
pflags  = false(1,nargs); n_pflags = 0;
switch length(sz),
    case 0, bands = [1,1];
    case 1, bands = sz * [1,1];
    otherwise, bands = sz(1:2);
end
for k = 2 : nargs,
    amin = 0; amax = 0;
    nm = lower( name{k} );
    if ~isempty(uplo),
        nm = [ uplo, '_', nm ]; %#ok
        name{k} = nm;
        orig{k} = [ uplo, '_', orig{k} ];
        uplo = '';
    end
    if nm(end) == '_',
        stypes{end+1} = nm; %#ok
        continue;
    end
    switch nm,
        case { 'binary', 'integer' },
            itype = nm;
            iflags(k) = true;
        case { 'upper', 'lower', 'skew', 'doubly' },
            uplo = nm;
            continue;
        case 'nonnegative',
            stypes{end+1} = nm; %#ok
            do_nneg = true;
        case 'complex',        
            do_comp = true;
            if do_semi, do_conj = true; end
            do_strx = true;
        case 'complex_if',
            amin = 1; amax = 1;
            if length(args{k}) == 1,
                arg = args{k}{1};
                if isa(arg,'logical') && arg || ~isreal(arg),
                    do_comp = true;
                    if do_semi, do_conj = true; end
                    do_strx = true;
                end
            end
        case 'symmetric',
            sflags(k) = true;                 do_symm = true;
            do_strx = true;
        case 'symmetric_ut',
            sflags(k) = true;
            pflags(k) = true;
            do_strx = true;
        case 'hermitian',      
            sflags(k) = true; do_comp = true; do_symm = true; do_conj = true;
            do_strx = true;
        case 'hermitian_if',
            amin = 1; amax = 1;
            sflags(k) = true; do_symm = true;
            if length(args{k}) == 1,
                arg = args{k}{1};
                if isa(arg,'logical') && arg || ~isreal(arg),
                    do_comp = true; do_conj = true;
                    do_strx = true;
                end
            end
        case { 'skew_symmetric', 'skew-symmetric' },
            sflags(k) = true;                 do_symm = true;                 do_skew = true;
            do_strx = true;
        case { 'skew_hermitian', 'skew-hermitian' }
            sflags(k) = true; do_comp = true; do_symm = true; do_conj = true; do_skew = true;
            do_strx = true;
        case {'hankel','upper_hankel'},
            pflags(k) = true;
            sflags(k) = true;
            hflags(k) = true;
            do_strx = true;
        case 'toeplitz',
            hflags(k) = true;
            is_toep = true;
            do_strx = true;
        case { 'semidefinite', 'doubly_nonnegative' },
            do_semi = true; do_symm = true;
            if do_comp, do_conj = true; end
            stypes{end+1} = 'semidefinite'; %#ok
            if nm(1) == 'd', 
                stypes{end+1} = 'nonnegative'; %#ok
                do_nneg = true; do_dble = true; 
            end
            do_strx = true;
        case 'upper_bidiagonal',
            bflags(k) = true; sflags(k) = true;
            bands = min(bands,[0,1]);
            do_strx = true;
        case 'lower_bidiagonal',
            bflags(k) = true; sflags(k) = true;
            bands = min(bands,[1,0]);
            do_strx = true;
        case 'upper_hessenberg',
            bflags(k) = true; sflags(k) = true;
            bands = min(bands,[1,Inf]);
            do_strx = true;
        case 'lower_hessenberg',
            bflags(k) = true; sflags(k) = true;
            bands = min(bands,[Inf,1]);
            do_strx = true;
        case 'upper_triangular',
            bflags(k) = true; sflags(k) = true;
            bands = min(bands,[0,Inf]);
            do_strx = true;
        case 'lower_triangular',
            bflags(k) = true; sflags(k) = true;
            bands = min(bands,[Inf,0]);
            do_strx = true;
        case 'diagonal',
            bflags(k) = true;
            bands = min(bands,[0,0]);
            do_strx = true;
        case 'tridiagonal',
            bflags(k) = true;
            bands = min(bands,[1,1]);
            do_strx = true;
        case 'scaled_identity',
            bflags(k) = true; is_toep = true;
            bands = min(bands,[0,0]);
            do_strx = true;
        case 'banded',
            amin = 1; amax = 2;
            sflags(k) = length(args{k}) == 2 && ~isequal(args{k}{1},args{k}{2});
            bflags(k) = true;
            lb = args{k}{1};
            switch length(args{k}),
                case 0, cvx_throw( 'Not enough arguments for "banded()" structure.' );
                case 1, ub = args{k}{1};
                case 2, ub = args{k}{2};
                otherwise, cvx_throw( 'Too many arguments for "banded()" structure.' );
            end
            if ~isnumeric( lb ) || length( lb ) ~= 1 || lb < 0 || lb ~= floor( lb ),
                cvx_throw( 'Bandwidth arguments must be nonnegative integers.' );
            elseif ~isnumeric( ub ) || length( ub ) ~= 1 || ub < 0 || ub ~= floor( ub ),
                cvx_throw( 'Bandwidth arguments must be nonnegative integers.' );
            end
            bands = min(bands,[lb,ub]);
            do_strx = true;
        otherwise,
            pflags(k) = true;
            if ~exist( [ 'cvx_s_', nm ], 'file' ),
                cvx_throw( 'Undefined matrix structure type: %s\nTrying to declare multiple variables? Use the VARIABLES keyword instead.', orig{k} );
            end
    end
    n_bflags = n_bflags + bflags(k);
    n_iflags = n_iflags + iflags(k);
    n_sflags = n_sflags + sflags(k);
    n_hflags = n_hflags + hflags(k);
    n_pflags = n_pflags + pflags(k);
    if length( args{k} ) < amin,
        cvx_throw( 'Not enough arguments: %s', orig{k} );
    elseif length( args{k} ) > amax,
        cvx_throw( 'Too many arguments: %s', orig{k} );
    end
end
if ~isempty(uplo),
    cvx_throw( 'Invalid structure type: %s', orig{end} );
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify integer consistency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_iflags,
    if n_iflags > 1,
        cvx_throw( 'Multiple integer keywords used:\n   %s', sprintf(' %s', orig{iflags} ) );
    elseif do_comp,
        xflags = iflags | findstr( name, 'complex', 'hermitian', 'skew_hermitian' );
        cvx_throw( 'Integer variables must also be real:\n   %s', sprintf(' %s', orig{xflags}) );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quick exit for no structure %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~do_strx && nargout > 0
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify structure consistency %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xflags = 0;
if n_sflags > 1, xflags = xflags | sflags; end
if n_hflags > 1, xflags = xflags | hflags; end
if n_iflags > 1, xflags = xflags | iflags; end
if n_bflags > 1, xflags = xflags | bflags; end
if do_semi
    if do_skew
        xflags = xflags | sflags | findstr( name, 'semidefinite', 'doubly_nonnegative' );
    elseif do_nneg && ~do_dble,
        xflags = xflags | findstr( name, 'nonnegative', 'semidefinite' );
    end
end
if do_comp && ( n_iflags || do_nneg )
    xflags = xflags | iflags | findstr( name, 'complex', 'hermitian', 'nonnegative', 'doubly_nonnegative' );
end
if do_nneg && do_skew
    xflags = xflags | findstr( name, 'nonnegative', 'skew_symmetric' );
end
if any( xflags ),
    cvx_throw( 'These forms of structure may not be specified simultaneously:\n   %s', sprintf(' %s', orig{xflags} ) );
end
if do_symm && sz(1) ~= sz(2),
    xflags = sflags | findstr( name, 'semidefinite', 'doubly_nonnegative' );
    cvx_throw( 'This matrix structure requires square matrices:%s', sprintf(' %s', orig{xflags} ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the basis matrices for the remaining structure elements        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strs = {};
if is_toep || any( bands < sz(1:2) ),
    [ strs{end+1}, do_symm ] = cvx_s_banded( sz(1), sz(2), [ do_symm, is_toep ], bands(1), bands(2) );
    do_symm = do_symm(1);
end
if n_pflags,
    for k = find(pflags)',
        try
            [ strs{end+1}, do_symm ] = feval( [ 'cvx_s_', lower(name{k}) ], sz( 1 ), sz( 2 ), do_symm, args{k}{:} ); %#ok
        catch exc
            cvx_throw( 'Error constructing structure: %s\n   %s', orig{k}, exc.message );
        end
    end
end
if do_symm,
    [ strs{end+1}, do_symm ] = cvx_s_symmetric( sz(1), sz(2), do_symm ); %#ok
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If multiple structures have been requested (e.g., toeplitz and banded),  %
% combine them together by finding bases for their orthogonal complements, %
% concatenating, and taking the orthogonal complement of that. This should %
% be used much less frequently than before---if ever---now that we handle  %
% symmetry as a special case for improved performance.                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch length( strs ),
    case 0,
        if ~do_comp && ~do_skew && nargout ~= 0, return; end
        sz(end+1:2) = 1;
        nel = sz( 1 ) * sz( 2 );
        S = sparse( 1 : nel, 1 : nel, 1, nel, nel );
    case 1,
        S = strs{ 1 };
    otherwise,
        for k = 1 : length(strs),
            strs{k} = cvx_orthog_structure( strs{k} ); %#ok
        end
        S = cvx_orthog_structure( vertcat(strs{:}), true );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle complex, skew-symmetric, and Hermitian structures.                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if do_comp,
    S = [ S ; +1j * S ];
    S = S([1:end/2;end/2+1:end],:);
end
if do_skew || do_conj,
    r = (0:sz(1)-1)'; r = r(:,ones(1,sz(2)));
    c = 0:sz(2)-1; c = c(ones(1,sz(1)),:);
    ut = r < c; dg = r == c;
    if do_skew,
        S(:,ut) = - S(:,ut);
        S(:,dg) = 0;
    end
    if do_conj,
        S(:,ut) = conj(S(:,ut));
        S(:,dg) = real(S(:,dg));
    end
    S = S(any(S,2),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Report an error of the structure is empty                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty( S ),
    cvx_throw( 'Incompatible structure modifiers:%s', sprintf( ' %s', args.orig ) );
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replicate structure for N-D arrays                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length( sz ) > 2,
    S = cvx_replicate_structure( S, sz( 3 : end ) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the structure if called with no output arguments                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 0,
    [ii,jj,vv] = find( S );
    if isempty(sz), sz = [1,1]; end
    Z = reshape( full( sparse( jj, 1, ii .* vv, prod(sz), 1 ) ), sz );
    temp = sprintf( ',%d', sz );
    fprintf( 'Name: %s\nSize: %s\nStructure:', name{1}, temp(2:end) );
    if ~isempty( itype ),
        stypes = [ stypes, itype ];
    end
    for k = 2 : nargs,
        if ~any( strcmp(orig{k},stypes) )
            fprintf( ' %s', orig{k} );
        end
    end
    fprintf( '\n' );
    if ~isempty( stypes ),
        fprintf( 'Passthrough:' );
        for k = 1 : length(stypes),
            fprintf( ' %s', stypes{k} );
        end
    end
    fprintf( '\n\n' );
    fmt = get(0,'format');
    set(0,'format','rational');
    disp( Z );
    set(0,'format',fmt);
    clear S
end

function z = findstr( x, varargin )
q = varargin;
z = cellfun( @(z)any(strcmpi(z,q)), x );

% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for orig copyright information.
% The command 'cvx_where' will show where this file is located.
