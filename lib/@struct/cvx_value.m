function x = cvx_value( x )
x = cell2struct( cvx_value( struct2cell( x ) ), fieldnames( x ), 1 );

