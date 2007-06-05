function tf = cvx_use_sparse( sz, nz, isr )
if nargin == 1,
    ss = size( sz );
    if any( ss == 1 ) | length( ss ) > 2,
        tf = false;
        return
    end
    isr = isreal( sz );
    if issparse( sz ),
        nz = nzmax( sz );
    else
        nz = nnz( sz );
    end
    sz = ss;
elseif any( sz == 1 ) | length( sz ) > 2,
    tf = false;
    return
elseif nargin < 3,
    isr = true;
end
if isr,
    tf = 1 + ( 1 - 2 * sz( 1 ) ) * sz( 2 ) + 3 * nz < 0;
else
    tf = 1 + ( 1 - 4 * sz( 1 ) ) * sz( 2 ) + 5 * nz < 0;
end
