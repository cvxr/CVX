% PENALTY_COMP_CVX         Penalty function approximation
%                          (a figure is generated)
% Sec. 6.1.2, fig 6.2, Boyd & Vandenberghe "Convex Optimization"
% Original by Lieven Vandenberghe 
% Adapted for CVX Argyris Zymnis - 10/2005
%
% Comparison of the ell1, ell2, deadzone-linear and log-barrier
% penalty functions for the approximation problem:
%       minimize phi(A*x-b),
%
% where phi(x) is the penalty function
% Log-barrier will be implemented in the future version of CVX

dz=0.5;  % scaling for deadzone penalty

% pick random mxn A and m-vector b and scale so that |b_i| < 1, so 
% that x=0 is feasible for log-barrier
randn('seed',0);
m=100;  n=30;  
A=randn(m,n);
b=randn(m,1); 

% ell_1 approximation
% 
% minimize   ||Ax+b||_1

disp('ell-one approximation');
cvx_begin
    variable x1(n)
    minimize(norm(A*x1+b,1))
cvx_end


% ell_2 approximation
%
% minimize ||Ax+b||_2

disp('ell-2');
x2=-A\b;


% deadzone penalty approximation 
%
% minimize   1'*y
% subject to -y <= Ax+b <= y
%             y >= dz

disp('deadzone penalty');
cvx_begin
    variable xdz(n)
    minimize(sum(deadzone(A*xdz+b,dz)))
cvx_end


% log-barrier penalty approximation 
%
% minimize -sum log(1-(ai'*x+bi)^2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cvx_begin
%     variable xlb(n)
%     minimize(-sum(log(1-square(A*xlb+b))))
% cvx_end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('log-barrier')

% parameters for Newton Method & line search
alpha=.01; beta=.5;  

% minimize linfty norm to get starting point 
cc = [zeros(n,1); 1];  
AA = [A -ones(m,1); -A -ones(m,1)]; 
bb = [-b; b];
x = linprog(cc,AA,bb);
linf = x(n+1);
xlb = x(1:n);
A = A/(1.1*linf);
b = b/(1.1*linf);

for iters = 1:50

   yp = 1 - (A*xlb+b);  ym = (A*xlb+b) + 1;
   f = -sum(log(yp)) - sum(log(ym));  
   g = A'*(1./yp) - A'*(1./ym);
   H = A'*diag(1./(yp.^2) + 1./(ym.^2))*A;
   v = -H\g;
   fprime = g'*v;
   ntdecr = sqrt(-fprime);
   if (ntdecr < 1e-5), break; end;  

   t = 1;
   newx = xlb + t*v;
   while ((min(1-(A*newx +b)) < 0) | (min((A*newx +b)+1) < 0)) 
       t = beta*t;
       newx = xlb + t*v;
   end;
   newf = -sum(log(1 - (A*newx+b))) - sum(log(1+(A*newx+b)));
   while (newf > f + alpha*t*fprime)
       t = beta*t;
       newx = xlb + t*v;
       newf = -sum(log(1-(A*newx+b))) - sum(log(1+(A*newx+b)));
   end;
   xlb = xlb+t*v; 
end


% Plot histogram of residuals

ss = max(abs([A*x1+b; A*x2+b; A*xdz+b;  A*xlb+b]));
tt = -ceil(ss):0.05:ceil(ss);  % sets center for each bin
[N1,hist1] = hist(A*x1+b,tt);
[N2,hist2] = hist(A*x2+b,tt);
[N3,hist3] = hist(A*xdz+b,tt);
[N4,hist4] = hist(A*xlb+b,tt); 


range_max=2.0;  rr=-range_max:1e-2:range_max;

figure(1), clf, hold off
subplot(4,1,1), 
bar(hist1,N1);  
hold on
plot(rr, abs(rr)*40/3, '-');
ylabel('p=1')
axis([-range_max range_max 0 40]);
hold off

subplot(4,1,2), 
bar(hist2,N2); 
hold on;
plot(rr,2*rr.^2), 
ylabel('p=2') 
axis([-range_max range_max 0 11]);
hold off

subplot(4,1,3), 
bar(hist3,N3); 
hold on
plot(rr,30/3*max([zeros(1,length(rr));abs(rr)-dz]))
ylabel('Deadzone')
axis([-range_max range_max 0 25]);
hold off

subplot(4,1,4), 
bar(hist4,N4); 
rr_lb=linspace(-1+(1e-6),1-(1e-6),600);
hold on
plot(rr_lb, -3*log(1-rr_lb.^2),rr,2*rr.^2,'--')
axis([-range_max range_max 0 11]);
ylabel('Log barrier'), 
xlabel('r')
hold off
