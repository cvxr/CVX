function prob  = cher(A,b,Sigma);

% calculates upper bound that Gaussian random vector N(0,Sigma) 
% satisfies Ax <= b
%
% minimize  b'*u + (1/2) u'*A*Sigma*A'*u 
% s.t       u >= 0
%
% solve as an SDP
% 
% minimize b'*u + t/2
% s.t.     [ I       R*A'*u ]
%          [ u'*A*R'  t     ] > = 0
%          u >= 0
%


m = size(A,1);
novars = 1+m;
R = chol(Sigma);

% minimize b'*u + t/2
c = [b; 1/2];

szs = [3;ones(m,1)];
F = zeros(sum(szs.^2), novars+1);


% block 1
blck = [eye(2) zeros(2,1); zeros(1,2) 0 ];
F(1:9,1) = blck(:);
for i=1:m
   blk = zeros(3,3);
   blk(1:2,3) = R*A(i,:)';
   blk(3,1:2) = A(i,:)*R';
   F(1:9, i+1) = blk(:);
end;
blk = [0 0 0; 0 0 0; 0 0 1];
F(1:9, m+2) = blk(:);

% scalar blocks
F(9+[1:m], 2:(m+1)) = eye(m);


% initial point
u = ones(m,1);
t = norm(R*A'*u)^2 + 1;
x0 = [u;t];

[x,Z,z,ul,time] = bigM(F,szs,c,x0,1e5,10,1e-8,1e-4,0,100);
prob = exp(ul(1));

