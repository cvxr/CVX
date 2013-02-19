function S = cvx_create_structure( sz, varargin )

%CVX_CREATE_STRUCTURE Construct a basis for a matrix structure.

error( nargchk( 1, Inf, nargin ) ); %#ok

[ temp, sz ] = cvx_check_dimlist( sz, false );
if ~temp,
    error( 'First argument must be a non-empty dimension list.' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scan the structure strings for symmetry and complex cases. We now handle %
% these cases specially for improved performance.                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do_skew = false;
do_conj = false;
do_comp = false;
do_symm = false;
is_toep = false;
is_hank = false;
nstrs = nargin - 1;
oargs = varargin;
strs = {};
for k = 1 : nstrs,
    argt = varargin{k};
    if isstruct( argt ),
        strx = argt;
        argt = argt.name; 
    elseif any( argt == '(' ),
        if argt(end) ~= ')',
            error( 'CVX:InvalidStructure', 'Invalid matrix structure: %s', arg );
        end
        xt = find( argt == '(' );
        strx.name = argt(1:xt(1)-1);
        strx.args = evalin( 'caller', [ '{', argt(xt(1)+1:end-1), '}' ] );
        strx.full = argt;
        varargin{k} = strx;
        argt = strx.name;
    end
    switch lower( argt ),
        case 'complex',        do_comp = true; varargin{k} = [];
        case 'symmetric',      strs{end+1} = argt; do_symm = true; varargin{k} = []; %#ok
        case 'hermitian',      strs{end+1} = argt; do_symm = true; do_comp = true; do_conj = true; varargin{k} = []; %#ok
        case 'skew_symmetric', strs{end+1} = argt; do_symm = true; do_skew = true; varargin{k} = []; %#ok
        case 'skew_hermitian', strs{end+1} = argt; do_symm = true; do_skew = true; do_comp = true; do_conj = true; varargin{k} = []; %#ok
        case {'hankel','upper_hankel'}, strs{end+1} = argt; is_hank = true; %#ok
        case 'toeplitz', is_toep = true;
        case 'symmetric_ut', strs{end+1} = argt; %#ok
        case {'upper_bidiagonal','upper_triangular','upper_hessenberg',...
              'lower_bidiagonal','lower_triangular','lower_hessenberg'},
            strs{end+1} = argt; %#ok
        case 'banded',
            if length(strx.args) == 2 && ~isequal(strx.args{1},strx.args{2}),
                strs{end+1} = strx.full; %#ok
            end
        otherwise,
            name = [ 'cvx_s_', lower(argt) ];
            if ~exist( name, 'file' ),
                error( 'CVX:InvalidStructure', 'Undefined matrix structure type: %s', argt );
            end
    end
end
symm = length( strs );
if symm,
    if symm > 1,
        error( 'CVX:InvalidStructure', 'These forms of structure may not be specified simultaneously:\n   %s', sprintf(' %s',strs{:}) );
    elseif is_toep && is_hank,
        error( 'CVX:InvalidStructure', 'Cannot specify both Toeplitz and Hankel structure.' );
    elseif sz(1) ~= sz(2),
        error( 'CVX:InvalidStructure', 'Symmetric structure square matrices.' );
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the basis matrices for the remaining structure elements        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strs = {};
for k = 1 : nstrs,
    argt = varargin{k};
    if isempty( argt ),
        continue;
    elseif isstruct( argt ),
        try
            args = argt.args;
        catch %#ok
            args = {};
        end
        argt = argt.name;
    else
        args = {};
    end
    [ S, do_symm ] = feval( [ 'cvx_s_', lower(argt) ], sz( 1 ), sz( 2 ), do_symm, args{:});
    strs{ end + 1 } = S; %#ok
end
if do_symm,
    [ strs{ end + 1 }, do_symm ] = cvx_s_symmetric( sz(1), sz(2), do_symm ); %#ok
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
    for k = 1 : length(oargs),
        if isstruct(oargs{k}),
            oargs{k} = oargs{k}.full;
        end
        error( 'CVX:InvalidStructure', 'Incompatible structure modifiers:%s', sprintf( ' %s', oargs{:} ) );
    end
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
    Z = reshape( full( sparse( jj, 1, ii .* vv, prod(sz), 1 ) ), sz );
    temp = sprintf( ',%d', sz );
    fprintf( '\n(%s)', temp(2:end) );
    for k = 1 : length(oargs),
        argt = oargs{k};
        if isstruct(argt), argt = argt.full; end
        fprintf( ' %s', argt );
    end
    fprintf( '\n\n' );
    disp( Z );
    clear S
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
