function s = cvx_profile( flag )
global cvx___
if isempty( cvx___ ), 
    cvx_setpath( 1 ); 
end
s = cvx___.profile;
if nargin == 1,
    if length( flag ) ~= 1,
        error( 'Argument must be a numeric or logical scalar.' );
    end
    flag = double( flag ) ~= 0;
    if cvx___.profile ~= flag,
        cvx___.profile = flag;
        stat = profile('status');
        if flag & ~isempty( cvx___.problems ) & ~isequal( stat.ProfilerStatus, 'on' ),
            profile on
        elseif ~flag & isequal( stat.ProfilerStatus, 'on' ),
            profile off
        end
    end
end
    

