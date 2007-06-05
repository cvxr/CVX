function v = cvx_value( x )
warning( sprintf( 'CVX error: illegal use of a cvx problem object has been detected.\n   Please do not copy or manipulate the value of ''cvx_problem'' in any way.' ) );
v = [];

