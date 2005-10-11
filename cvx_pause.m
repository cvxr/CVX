function s = cvx_pause( flag )
global cvx___
if isempty( cvx___ ), cvx_setpath( 1 ); end
s = cvx___.pause;
if nargin == 1,
    if length( flag ) ~= 1,
        error( 'Argument must be a numeric or logical scalar.' );
    end
    cvx___.pause = double( flag ) ~= 0;
end
