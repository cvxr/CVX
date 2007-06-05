function x = cvx_zeros( s )
if cvx_use_sparse( s, 0, 1 ),
     x = sparse( s(1), s(2) );
else
     x = zeros( s );
end
