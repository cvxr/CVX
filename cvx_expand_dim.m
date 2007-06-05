function zx = cvx_expand_dim( z, dim, nx )
zdims = cell( 1, max( ndims(z), dim ) );
[ zdims{:} ] = deal( ':' );
zdims{dim} = ones( 1, nx );
zx = z( zdims{:} );

