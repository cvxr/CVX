function cvx_optpnt = nonnegative( sx ) %#ok

%NONNEGATIVE   The nonnegative orthant.
%   NONNEGATIVE(SX), where SX is a valid size vector, creates an array
%   of size SX and constrains each element to be nonnegative. Therefore,
%   given the declaration
%      variable x(sx)
%   the constraint
%      x == nonnegative(sx);
%   is equivalent to
%      x >= 0;
%   Obviously, the inequality form is simpler and is preferred in most
%   circumstances.
%
%   Disciplined convex programming information:
%       NONNEGATIVE is a cvx set specification. See the user guide for
%       details on how to use sets.

error( nargchk( 1, 1, nargin ) );

[ temp, sx ] = cvx_check_dimlist( sx, true );
if ~temp,
    error( 'Argument must be a non-empty dimension vector.' );
end
    
cvx_begin_set
    variables x( sx )
    if all( sx ~= 0 ),
        [ tx, dummy ] = find( cvx_basis( x ) ); %#ok
        newnonl( cvx_problem, 'nonnegative', tx(:) );
    end
cvx_end_set

% Copyright 2010 Michael C. Grant and Stephen P. Boyd. 
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
