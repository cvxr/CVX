function gdim = cvx_pushexp( ndim )

global cvx___
persistent expv geov
if isempty( expv ),
    expv = int8([3,3,3,22,15,15,15,16,16,16,17,17,17,22,22,17,17,17,21,21,21,22])';
    geov = int8([3,3,3,22,22,22,22,19,19,19,21,21,21,22,22,19,19,19,21,21,21,22])';
end
if cvx___.geometric,
    map = geov;
else
    map = expv;
end
gdim = cvx_pushvar( numel(ndim), map(cvx___.classes(ndim)) );
cvx___.logarithm( gdim, 1 ) = ndim(:);
cvx___.exponential( ndim, 1 ) = gdim(:);
