function s = cvx_profile( flag )

% CVX_PROFILE	CVX-specific profiler control.
%    This is a function used for internal CVX development to help determine 
%    performance limits within the CVX code itself, by turning off the profiler
%    when the solver is being called. End users will likely not find this
%    function to be useful.

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
  
% Copyright 2007 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
  

