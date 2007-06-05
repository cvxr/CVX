function y = cvx_value( x )
global cvx___
if cvx___.mversion < 7.1,
    sx = size( x );
    y = cell( sx );
    for k = 1 : prod( sx ),
        y{k} = cvx_value( x{k} );
    end
else
    y = cellfun( @cvx_value, x, 'UniformOutput', false );
end

