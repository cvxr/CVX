function varargout = cvx_get_dimlist( args )
if isempty( args ),
	sx = [1,1];
else
    if iscell( args ),
        sx = args{1};
    else
        sx = args;
    end
    nx = numel(sx);
    if ~( isnumeric(sx) && nx==size(sx,2) && isreal(sx) && all(sx>=0) && all(sx==floor(sx)) ),
        cvx_throw( 'Size argument must be a row vector of nonnegative integers.' );
    end
    switch nx,
    case 2,
    case 1,
        sx = [ sx ,1 ];
    otherwise,
        if sx(end) == 1,
            nx = max([2,find(sx~=1,1,'last')]);
            sx = sx(1:nx);
        end
    end
end
if iscell( args ),
    args{1} = sx;
else
    args = { sx };
end
if nargout > length(args), 
    args{nargout} = []; 
end
varargout = args;
