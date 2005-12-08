%
% Maximum entropy reconstruction example
%
% distribution of n points, uniform between -1 and 1 
%

n = 100;  a = linspace(-1,1,n)';

%moment constraints Ax <= b
A = sparse([a'; -a';                      
   (a.^2)';  -(a.^2)';
   (3*a.^3-2*a)'; -(3*a.^3-2*a)';  
   (a<0)';  -(a<0)'; 
   -eye(n)]);
b = [0.1; 0.1; 
      0.6; -0.5;
     -0.2; 0.3;
     0.4; -0.3; 
     zeros(n,1)];


% eliminate last variable
AA = A(:,1:n-1) - A(:,n)*ones(1,n-1);
bb = b - A(:,n);
bb(length(bb)) = 1;


% maximize and minimize cumsum  up to t, t=1,...,100
Ubnds = zeros(1,n);  Ubnds(n) = 1;
Lbnds = zeros(1,n);  Lbnds(n) = 1;


for t=1:(n-1)

   t

   cc = [ones(t,1); zeros(n-t-1,1)];

   %[z,x,status] = mpc(AA', cc, bb);
   x = linprog(cc,AA,bb);
   Lbnds(t) = cc'*x;

   x = linprog(-cc,AA,bb);
   Ubnds(t) = cc'*x;

end



% 
% max entropy distribution
%
% min  sum_i xi*log xi
% s.t Ax <= b;
%     1'*x = 1
%     x >= 0
%

%
% feasibility problem for initial point
%
 
cc = [zeros(n-1,1);1];
AA = [AA -ones(size(AA,1),1)];
%[z,x,status] = mpc(AA', -cc, bb);
xx = linprog(cc,AA,bb);
x = [xx(1:n-1); 1-sum(xx(1:n-1))];

slack = min(b-A*x);
if (slack < 0) keyboard; end;


MU = 5;
ALPHA = 0.1;
BETA = 0.5;
TOL = 1e-4;
MAXITERS = 50;

m = size(A,1);

t = 1;

while ((size(A,1))/t > TOL)

   for i=1:MAXITERS

       d = b-A*x;
       val = t*x'*log(x) - sum(log(d));
       g = t*(1+log(x)) + A'*(1./d);
       H = t*diag(1./x) + A'*diag(1./(d.^2))*A;

       sol = [H ones(n,1); ones(1,n) 0] \ [-g; 0];      
       v = sol(1:n);
       fprime = g'*v;
       if (abs(fprime) < 1e-5) break; end;

       ntdecr = sqrt(-fprime)

       s = 1;
       newx = x + s*v;
       while ((min(newx) < 0)  | (min(b-A*newx) < 0)) 
           s = s*BETA;
           newx = x + s*v;
       end;
       newx = x + s*v;
       newval = t*newx'*log(newx) - sum(log(b-A*newx));
       while (newval > val + s*ALPHA*fprime)
           s = s*BETA;
           newx = x + s*v;
           newval = t*newx'*log(newx) - sum(log(b-A*newx));
        end;

        x =  x + s*v;

   end;

   t = MU*t;

end;

figure(1)
xent = x; 
stairs(a,cumsum(xent));
grid on;
hold on
d = stairs(a, Lbnds,'r-');  set(d,'Color',[0 0.5 0]);
d = stairs(a, Ubnds,'r-');  set(d,'Color',[0 0.5 0]);  
d = plot([-1,-1], [Lbnds(1), Ubnds(1)],'r-'); 
set(d,'Color',[0 0.5 0]);  
axis([-1.1 1.1 -0.1 1.1]);
%xlabel('x');
ylabel('y');
hold off

%axis([-1 1 0 1.1*max(xent)]);
xlabel('x');
ylabel('y');
hold off
print -deps maxentcumsum.eps


figure(2);
stairs(a,xent);
xlabel('x');
ylabel('y');
print -deps maxentdistr.eps
