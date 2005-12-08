% BASISPURSUIT   Basis pursuit using Gabor functions
%                (figures are generated)
% Sec. 6.5.4, fig 6.21, 6.22 and 6.23 
% Boyd & Vandenberghe "Convex Optimization"
% Original by Lieven Vandenberghe 
% Adapted for CVX by Argyris Zymnis - 11/27/2005
% 
% Here we find a sparse basis for a signal y out of 
% a set of Gabor functions. We do this by solving
%       minimize  ||A*x-y||_2 + ||x||_1
%
% where the columns of A are sampled Gabor functions.
% We then fix the sparsity pattern obtained and solve
%       minimize  ||A*x-y||_2
%
% NOTE: The file takes a while to run

clear

% Problem parameters
sigma = 0.05; % Size of Gaussian function
T = 0.002; % Sample time
Thr = 0.001; % Basis signal threshold
kmax = 30; % Number of signals are 2*kmax+1
w0 = 5;    % Base frequency (w0 * kmax should be 150 for good results)

% Build dictionary matrix
disp('Building dictionary matrix...');
A = sparse(401,401*61);
p = 1;
t = (0:T:1)';
for s = 0:T:1
  x = exp(-(t-s).^2/(sigma^2)); x(find(x < Thr)) = 0;
  A(find(x),p) = x(find(x)); % Gaussian
  p = p+1;
  for k=1:kmax
    c = x.*sin(w0*t*k);
    A(find(c),p) = c(find(c)); % Sinus basis
    p = p + 1;
    c = x.*cos(w0*t*k);
    A(find(c),p) = c(find(c)); % Cosine basis
    p = p + 1;
  end;
end;
disp('done.');

% Construct example signal
a = 0.5*sin(t*11)+1; 
theta = sin(5*t)*30;
b = a.*sin(theta);

% Solve the Basis Pursuit problem
disp('Solving Basis Pursuit problem...');
tic
cvx_begin
    variable x(30561)
    minimize(square_pos(norm(A*x-b,2))+norm(x,1))
cvx_end
disp('done');
toc

% Reoptimize problem over nonzero coefficients
p = find(abs(x) > 1e-5);
A2 = A(:,p);
x2 = A2 \ b;

% Constants
M = 61; % Number of different Basis signals
sk = 250; % Index of s = 0.5 

% Plot example basis functions;  
%if (0) % to do this, re-run basispursuit.m to create A
figure(1); clf;
subplot(3,1,1); plot(t,A(:,M*sk+1)); axis([0 1 -1 1]);
title('Basis function 1');
subplot(3,1,2); plot(t,A(:,M*sk+31)); axis([0 1 -1 1]);
title('Basis function 2');
subplot(3,1,3); plot(t,A(:,M*sk+61)); axis([0 1 -1 1]);
title('Basis function 3');
%print -deps bp-dict_helv.eps

% Plot reconstructed signal
figure(2); clf;
subplot(2,1,1);
plot(t,A2*x2,'--',t,b,'-'); axis([0 1 -1.5 1.5]);
xlabel('t'); ylabel('y_{hat} and y');
title('Original and Reconstructed signals')
subplot(2,1,2);
plot(t,A2*x2-b); axis([0 1 -0.06 0.06]);
title('Reconstruction error')
xlabel('t'); ylabel('y - y_{hat}');
%print -deps bp-approx_helv.eps

% Plot frequency plot
figure(3); clf;

subplot(2,1,1); 
plot(t,b); xlabel('t'); ylabel('y'); axis([0 1 -1.5 1.5]);
title('Original Signal')
subplot(2,1,2); 
plot(t,150*abs(cos(w0*t)),'--'); 
hold on;
for k = 1:length(t);
  if(abs(x((k-1)*M+1)) > 1e-5), plot(t(k),0,'o'); end;
  for j = 2:2:kmax*2
    if((abs(x((k-1)*M+j)) > 1e-5) | (abs(x((k-1)*M+j+1)) > 1e-5)),
      plot(t(k),w0*j/2,'o');
    end;
  end;
end;
xlabel('t'); ylabel('w');
title('Instantaneous frequency')
hold off;
