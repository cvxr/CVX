function gdim = cvx_pushexp( ndim, nonew )

global cvx___
persistent expv geov
if isempty( expv ),
    expv = int8([3,3,3,4,15,15,15,16,16,16,17,17,17,22,22,17,17,17,17,17,17,22])';
    geov = int8([3,3,3,4,22,22,22,19,19,19,20,20,20,22,22,20,20,20,20,20,20,22])';
end
if cvx___.problems(end).gp,
    map = geov;
else
    map = expv;
end
cx = cvx___.classes(ndim);
cy = map(cx);
if nargin < 2 || ~nonew,
    gdim = cvx_newvar( numel(ndim), cy );
    cvx___.logarithm( gdim, 1 ) = ndim(:);
    cvx___.exponential( ndim, 1 ) = gdim(:);
else
    gdim = cvx___.exponential( ndim );
    cvx___.classes( gdim ) = cy;
end

