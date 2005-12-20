function [prob,P,q,r,X,lambda] = cheb(A,b,Sigma);

% calculates lower bound on probability than random vector x in R2
% with mean zero and covariance Sigma satisfies Ax <= b
%
% Sigma must be positive definite
%
% output arguments:
% - prob: lower bound on probability
% - P,q,r: x'*P*x + 2*q'*x + r is a quadratic function
%   that majorizes the 0-1 indicator function of the complement
%   of the polyhedron, 
% - X, lambda:  a discrete distribution with mean zero, covariance 
%   Sigma and Prob(X not in C)  >= 1-prob

%
% maximize  1 - Tr Sigma*P - r 
% s.t.      [ P                 q-taui*ai/2     ] 
%           [ (q-taui*ai/2)'    r-1 + taui*b(i) ]  >= 0 i=1,...m
%           taui >= 0
%           [ P q  ]
%           [ q' r ] >= 0
%
% variables P in Sn, q in Rn, r in R
%

m = size(A,1);
novars = 3+2+1+m;  % P,q,r,tau

% minimize Tr Sigma*P + r
c = [Sigma(1,1); 2*Sigma(2,1); Sigma(2,2); 0; 0; 1; zeros(m,1)]; 

szs = [3*ones(m,1); ones(m,1); 3];
F = zeros(sum(szs.^2), novars+1);

% blocks 1...m;
for i=1:m
   blk = [0 0 0; 0 0 0; 0 0 -1];
   F((i-1)*9+[1:9],1) = blk(:); 
   blk = [1 0 0 ; 0 0 0; 0 0 0];  % coeff of P(1,1)
   F((i-1)*9+[1:9],2) = blk(:);
   blk = [0 1 0 ; 1 0 0; 0 0 0];  % coeff of P(2,1)
   F((i-1)*9+[1:9],3) = blk(:);
   blk = [0 0 0 ; 0 1 0; 0 0 0];  % coeff of P(2,2)
   F((i-1)*9+[1:9],4) = blk(:);
   blk = [0 0 1 ; 0 0 0; 1 0 0];  % coeff of q(1)
   F((i-1)*9+[1:9],5) = blk(:);
   blk = [0 0 0 ; 0 0 1; 0 1 0];  % coeff of q(2)
   F((i-1)*9+[1:9],6) = blk(:);
   blk = [0 0 0 ; 0 0 0; 0 0 1];  % coeff of r
   F((i-1)*9+[1:9],7) = blk(:);
   blk = (1/2)*[0 0 -A(i,1) ; 
                0 0 -A(i,2); 
               -A(i,1) -A(i,2) 2*b(i)];  % coeff of taui
   F((i-1)*9+[1:9],7+i) = blk(:);
end;

% taui >= 0, i=1,...m
F(m*9+[1:m], 7+[1:m]) = eye(m);

% last block
blk = [1 0 0 ; 0 0 0; 0 0 0];  % coeff of P(1,1)
F(10*m+[1:9],2) = blk(:);
blk = [0 1 0 ; 1 0 0; 0 0 0];  % coeff of P(2,1)
F(10*m+[1:9],3) = blk(:);
blk = [0 0 0 ; 0 1 0; 0 0 0];  % coeff of P(2,2)
F(10*m+[1:9],4) = blk(:);
blk = [0 0 1 ; 0 0 0; 1 0 0];  % coeff of q(1)
F(10*m+[1:9],5) = blk(:);
blk = [0 0 0 ; 0 0 1; 0 1 0];  % coeff of q(2)
F(10*m+[1:9],6) = blk(:);
blk = [0 0 0 ; 0 0 0; 0 0 1];  % coeff of r
F(10*m+[1:9],7) = blk(:);

% initial point
tau = ones(m,1);
t = 0;
for i=1:m
    blk = (1/2)*[0 0 tau(i)*A(i,1);  
                 0 0 tau(i)*A(i,2);  
                 tau(i)*A(i,1)   tau(i)*A(i,2)  2*(1-tau(i)*b(i)) ]; 
    t = max(t , max(eig(blk)));
end;
t = t+1;
P = t*eye(2);  q=[0;0];  r=t;
x0 = [P(1,1); P(1,2); P(2,2); q; r; tau];

[x,Z,z,ul,time] = bigM(F,szs,c,x0,1e4,10,1e-8,1e-6,0,100);
P = [x(1) x(2); x(2) x(3)];
q = [x(4); x(5)];
r = x(6);
prob = 1-ul(1);

X = [];
lambda = [];
for i=1:m
   Zi = reshape(Z((i-1)*9+[1:9]),3,3);
   if (abs(Zi(3,3)) > 1e-4)
      lambda = [lambda; Zi(3,3)];
      X = [X Zi(1:2,3)/Zi(3,3)];
   end;
end;
mu = 1-sum(lambda);
if (mu>1e-5)
   w = (-X*lambda)/mu;
   W = (Sigma - X*diag(lambda)*X')/mu;
   [v,d] = eig(W-w*w');
   d = diag(d);
   s = sum(d>1e-5);
   if (d(1) > 1e-5)
      X = [X w+sqrt(s)*sqrt(d(1))*v(:,1) ...
            w-sqrt(s)*sqrt(d(1))*v(:,1)];
      lambda = [lambda; mu/(2*s); mu/(2*s)];
   elseif (d(2) > 1e-5)
      X = [X w+sqrt(s)*sqrt(d(2))*v(:,2) ...
            w-sqrt(s)*sqrt(d(2))*v(:,2)];
      lambda = [lambda; mu/(2*s); mu/(2*s)];
   else
      X = [X w];
      lambda = [lambda; mu];
   end;
end;
