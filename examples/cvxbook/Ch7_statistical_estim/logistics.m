% logistics regresssion example

randn('state',0);
rand('state',0);

a =  1;
b = -5 ;
m= 100;

u = 10*rand(m,1);  
y = (rand(m,1) < exp(a*u+b)./(1+exp(a*u+b)));
plot(u,y,'o')
axis([-1,11,-0.1, 1.1]);

% max likelihood estimate
%
% minimize  -(sum_(y_i=1) ui)*a - b + sum log (1+exp(a*ui+b)
%

ind1 = find(y==1);
ind2 = find(y==0);
c = [-sum(u(ind1)); -length(ind1)];
A = [u ones(m,1)];

x = [1;-5];

ALPHA = 0.1;
BETA = 0.5;
for k=1:50

   val = c'*x + sum(log(1+exp(A*x)));
   w = exp(A*x)./(1+exp(A*x));
   g = c + A'*w;
   H = A'*diag(w./(1+exp(A*x)))*A;
   v = -H\g;

   fprime = g'*v
   ntdecr = sqrt(-fprime)

   if (ntdecr < 1e-8) break; end;
   s = 1;
   newx = x + s*v; 
   newval = c'*newx + sum(log(1+exp(A*newx)));
   while (newval > val + ALPHA*fprime*s) 
      s = BETA*s;
      newx = x + s*v; 
      newval = c'*newx + sum(log(1+exp(A*newx)));
   end;
   x = x+s*v;
end;

aml = x(1);  bml = x(2);
us = linspace(-1,11,1000)';
ps = exp(aml*us + bml)./(1+exp(aml*us+bml));

dots = plot(us,ps,'-', u(ind1),y(ind1),'o',...
     u(ind2),y(ind2),'o');

axis([-1, 11,-0.1,1.1]);
xlabel('x');
ylabel('y');
print -deps logistics.eps
