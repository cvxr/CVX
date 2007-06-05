function errs = cvx_lasterr()
errs = lasterr;
if strncmp( 'Error using ==>', errs, 15 ),
    errs = errs(min(find(errs==10))+1:end);
    while errs(end) == 10, errs(end) = []; end
end

