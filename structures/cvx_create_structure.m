function S = cvx_create_structure( sz, varargin )
%CVX_CREATE_STRUCTURE Construct a basis for a matrix structure.
error( nargchk( 1, Inf, nargin ) );

[ temp, sz ] = cvx_check_dimlist( sz, false );
if ~temp,
    error( 'First argument must be a non-empty dimension list.' );
elseif nargin <= 1,
    nel = prod( sz );
    S = sparse( 1 : nel, 1 : nel, 1, nel, nel );
    return
end

nstrs = nargin - 1;
strs  = cell( 1, nstrs );
rstr  = true;
for k = 1 : nstrs,
    argt = varargin{k};
    if isstruct( argt ),
        try
            args = argt.args;
        catch
            args = {};
        end
        argt = argt.name;
    else
        args = {};
    end
    name = [ 'cvx_s_', argt ];
    try
        if iscell( args ),
            S = feval( name, sz( 1 ), sz( 2 ), args{:} );
        else
            S = feval( name, sz( 1 ), sz( 2 ), args );
        end
    catch
        temp = exist( name, 'file' );
        if temp ~= 2 && temp ~= 3,
            error( 'Undefined matrix structure type: %s', argt );
        else
            error( lasterror );
        end
    end
    strs{ k } = S;
    if ~isreal( S ),
        rstr = false;
    end
end

if length( strs ) ~= 1,
    
    %
    % Each matrix in strs{:} contains the basis of a given matrix 
    % structure, arranged in rows. Here we intersect those bases by
    % finding bases for their orthogonal complements, combining them,
    % and taking the orthogonal complement of that. We then 'clean up'
    % the result by converting the result to reduced row echelon form.
    %
    
    n = size(strs{1},2)*(2-rstr);
    for k = 1 : nstrs,
        A = strs{k};
        if ~rstr,
            if isreal(A),
                A = cvx_blkdiag(A,A);
            else
                A = [real(A),imag(A)];
            end
            A = A(:,[1:n/2;n/2+1:n]);
        end
        strs{k} = cvx_orthog_structure(A);
    end
    NS = cvx_cleanup_structure( vertcat(strs{:}) );
    rr = size(NS,1);
    if rr == n,
        S = sparse(0,n);
    else
        S = cvx_orthog_structure(NS);
    end
    if ~rstr,
        S = S( :, 1 : 2: end ) + sqrt(-1) * S( :, 2 : 2 : end );
    end

end

if length( sz ) > 2,
    S = cvx_replicate_structure( S, sz( 3 : end ) );
end

% Copyright 2012 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
