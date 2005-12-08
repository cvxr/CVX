% POLY4_DISCR:  Polynomial discrimination 
%               (a figure is generated)
% Sec. 8.6.2, Boyd & Vandenberghe "Convex Optimization"
% Original by Lieven Vandenberghe 
% Adapted for CVX by Joelle Skaf - 10/23/05 

% The goal is to find the polynomial of degree 4 on IR^n that separates 
% two sets of points {x_1,...,x_N} and {y_1,...,y_N}. We are trying to find
% the coefficients of an order-4-polynomial P(x) that would satisfy: 
%           minimize    t
%               s.t.    P(x_i) <= t  for i = 1,...,N
%                       P(y_i) >= t   for i = 1,...,M

% Data generation 
rand('state',0);  randn('state',0);
N=100; 
M=120;
X1 = -1+2*rand(2,N);  X1 = X1*diag(0.9*rand(1,N)./sqrt(sum(X1.^2)));
ind = find(sqrt(sum((X1-[1.1*ones(1,N);zeros(1,N)]).^2)) < 0.9);
Y1 = X1(:,ind);
ind = find(sqrt(sum((X1-[1.1*ones(1,N);zeros(1,N)]).^2)) > 1);
X = X1(:,ind);
Y = -1+2*rand(2,M);  Y = Y*diag((1.1+rand(1,M))./sqrt(sum(Y.^2)));
Y = [Y Y1];

N = size(X,2);
M = size(Y,2);
% "Vandermonde"-like matrix 
monX = [ ones(1,N); X(1,:); X(2,:); X(1,:).^2;  X(1,:).*X(2,:); ...
         X(2,:).^2;  X(1,:).^3;  X(2,:).^2.*X(2,:); ...
         X(1,:).*X(2,:).^2;  X(2,:).^3;  X(1,:).^4; ...
         X(1,:).^3.*X(2,:); X(1,:).^2.*X(2,:).^2; X(1,:).*X(2,:).^3; ...
         X(2,:).^4 ];
monY = [ones(1,M); Y(1,:); Y(2,:); Y(1,:).^2;  Y(1,:).*Y(2,:); ...
        Y(2,:).^2;  Y(1,:).^3;  Y(2,:).^2.*Y(2,:); ...
        Y(1,:).*Y(2,:).^2;  Y(2,:).^3;  Y(1,:).^4; ...
        Y(1,:).^3.*Y(2,:); Y(1,:).^2.*Y(2,:).^2; Y(1,:).*Y(2,:).^3; ...
        Y(2,:).^4 ];

[m1,m2] = size(monX); 

% Solution via CVX
fprintf(1,'Finding the optimal polynomial of order 4 that separates the 2 classes...');

cvx_begin
    variables a(m1) t(1) 
    minimize ( t ) 
    a'*monX <= t
    a'*monY >= -t
cvx_end

fprintf(1,'Done! \n');

% Displaying results 
nopts = 2000;
angles = linspace(0,2*pi,nopts);
cont = zeros(2,nopts);
for i=1:nopts
   v = [cos(angles(i)); sin(angles(i))];
   l = 0;  u = 1;
   while (u-l > 1e-3)
      s = (u+l)/2;
      x = s*v;
      if (a'*[1; x(1); x(2); x(1)^2; x(1)*x(2); x(2)^2; x(1)^3; ...
             x(1)^2*x(2); x(1)*x(2)^2; x(2)^3; x(1)^4; x(1)^3*x(2); ...
             x(1)^2*x(2)^2; x(1)*x(2)^3; x(2)^4]  > 0), u = s;
      else, l=s; 
      end;
   end;
   cont(:,i) = s*v;
end;
  
graph = plot(X(1,:),X(2,:),'o', Y(1,:), Y(2,:),'o', cont(1,:), cont(2,:), '-');
set(graph(2),'MarkerFaceColor',[0 0.5 0]);
title('Optimal order-4 polynomial that separates the 2 classes')
% print -deps min-deg-discr.eps

%%%% Dual infeasible ????? 