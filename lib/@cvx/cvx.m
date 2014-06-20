function w = cvx( v, b )

%CVX   The CVX disciplined convex programming system.
%   CVX is a modeling framework for building, constructing, and solving
%   disciplined convex programs. CVX has online help for many functions and
%   operations, divided into several subsection:
%      cvx/commands   - Top-level commands to control CVX
%      cvx/keywords   - Keywords for declaring variables and objectives
%      cvx/builtins   - Built-in operators and functions supported in CVX models
%      cvx/functions  - Additional functions added by CVX
%      cvx/sets       - Definitions of common convex sets
%      cvx/structures - Matrix structure definitions
%   CVX also provides an extensive user guide in PDF format, which is found
%   in its top directory. This directory can be found by typing
%      cvx_where
%   at the command prompt.

global cvx___
switch nargin,
    case 2,
        w = cvx___.obj;
        switch numel( v ),
            case 2, 
                w.size_ = v;
            case 1, 
                w.size_ = [ v, 1 ];
            case 0,
                % w.size_ = [ 0, 1 ];
            otherwise,
                if v(end) == 1,
                    v = v(1:find(v>1,1,'last'));
                    v(1,end+1:2) = 1; 
                end
                w.size_ = v;
        end
        if isempty( b )
            w.basis_ = sparse( 1, prod( v ) );
        elseif issparse( b )
            w.basis_ = b;
        else
            w.basis_ = sparse( b );
        end
    case 1,
        if isnumeric( v ),
            w = cvx___.obj;
            w.size_ = size( v );
            if issparse( v ),
                w.basis_ = reshape( v, 1, prod(w.size_) );
            else
                w.basis_ = sparse( v(:)' );
            end
        elseif isa( v, 'cvx' ),
            w = v;
            return
        else
            w = class( v, 'cvx' );
        end
    case 0,
        w = cvx___.obj;
end
w.id_ = cvx___.id + 1;
cvx___.id = w.id_;
    
% Copyright 2005-2014 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
