function cvx_cd( subdir )

cd(cvx_where);
if nargin ~= 0,
    cd(subdir);
end

