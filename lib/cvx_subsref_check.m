function s = cvx_subsref_check( nin, nin_norm, S )

if nin < nin_norm,
    s = 'Not enough input arguments';
elseif nin == nin_norm,
    if ~isa( S, 'struct' ),
        s = 'Subscript argument to SUBSREF and SUBSASGN must be a structure.';
    else
        Sf = fieldnames( S );
        if length( Sf ) ~= 2,
            s = 'Subscript argument to SUBSREF and SUBSASGN must have two fields.';
        elseif ~isempty( setdiff( { 'type', 'subs' }, Sf ) ),
            s = 'Subscript argument to SUBSREF and SUBSASGN must have two fields whose names are "type" and "subs".';
        elseif isempty( S ),
            s = 'Subscript argument to SUBSREF and SUBSASGN must not be empty.';
        else
            s = '';
        end
    end
else
    s = '';
end

% Copyright 2010 Michael C. Grant and Stephen P. Boyd.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
