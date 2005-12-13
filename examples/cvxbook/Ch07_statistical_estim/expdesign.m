% Experiment design examples

m = 10;
angles1 = linspace(3*pi/4,pi,m);
angles2 = linspace(0,-pi/2,m);

V = [3.0*[cos(angles1); sin(angles1)], ...
     1.5*[cos(angles2); sin(angles2)]];
p = size(V,2);

noangles = 5000;

% plot sensor positions
figure(1)
plot(V(1,:),V(2,:),'ob',0,0,'k*');
axis([-4.5 4.5 -4.5 4.5])
set(gca,'XTick',[]);
set(gca,'YTick',[]);
axis off
%print -deps sensors.eps


% D-design


% max  log det V*diag(x)*V
% s.t.  sum(x)=1,  x >=0

MU = 10;
ALPHA = 0.01;
BETA = 0.5;
TOL = 1e-4;

t = 1;
x = (1/p)*ones(p,1);

while (p/t>TOL)

   for i=1:100

      % minimize t*logdet(V*diag(x)*V')^{-1} - sum log xi
      %   s.t.   1'*x = 1

      X = V*diag(x)*V';  R = chol(X);   % X = R'*R
      Xinv = inv(R)*inv(R');
      val = -2*t*sum(log(diag(R))) - sum(log(x));
      S = V'*Xinv*V;
      g = -t*diag(S) - 1./x;
      H = t*S.^2 + diag(1./(x.^2));
      sol = [H ones(p,1); ones(1,p) 0] \ [-g;0];
      v = sol(1:p);
      fprime = v'*g;
      ntdecr = sqrt(max(-fprime,0));
      if (ntdecr < 1e-4), break; end;
      s = 1;
      newx = x+s*v;  newX = V*diag(newx)*V';
      [newR, flag] = chol(newX);
      while ((flag) | min(newx) < 0),
          s = BETA*s;
          newx = x+s*v;  newX = V*diag(newx)*V';
          [newR, flag] = chol(newX);
      end;
      newx = x+s*v;  newX = V*diag(newx)*V';
      [newR, flag] = chol(newX);
      newval = -2*t*sum(log(diag(newR))) - sum(log(newx));
      while (newval > val + s*ALPHA*fprime)
          s = BETA*s;
          newx = x+s*v;  newX = V*diag(newx)*V';
          [newR, flag] = chol(newX);
          newval = -2*t*sum(log(diag(newR))) - sum(log(newx));
      end;

      x = x + s*v;

   end;

   t = MU*t;

end;


figure(1)
% draw ellipsoid v'*W*v <= 2
W = inv(V*diag(x)*V');
angles = linspace(0,2*pi,noangles);
R = chol(W);  % W = R'*R
ellipsoid = sqrt(2)*(R\[cos(angles); sin(angles)]);
d = plot(ellipsoid(1,:), ellipsoid(2,:), '--', 0,0,'+');
set(d, 'Color', [0 0.5 0]);
set(d(2),'MarkerFaceColor',[0 0.5 0]);
hold on;

dot=plot(V(1,:),V(2,:),'o');
ind = find(x > 0.001);
dots = plot(V(1,ind),V(2,ind),'o');
set(dots,'MarkerFaceColor','blue');
disp('Nonzero lambda values for D design:');
for i=1:length(ind)
   text(V(1,ind(i)),V(2,ind(i)), ['l',int2str(ind(i))]);
   disp(['lambda(',int2str(ind(i)),') = ', num2str(x(ind(i)))]);
end;

%axis([-4.5 4.5 -4.5 4.5])
axis([-5 5 -5 5])
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
hold off
axis off
%print -deps Ddesign.eps
xD = x;

keyboard


% A-design

%
% minimize Trace (sum_i xi*vi*vi')^{-1}
% s.t.     x >= 0, 1'*x = 1
%


MU = 10;
ALPHA = 0.1;
BETA = 0.5;
TOL = 1e-4;

t = 1;
x = (1/p)*ones(p,1);

while (p/t>TOL)

   for i=1:100

      % minimize t*trace(V*diag(x)*V')^{-1} - sum log xi
      %   s.t.   1'*x = 1

      X = V*diag(x)*V';  R = chol(X);   % X = R'*R
      Xinv = inv(R)*inv(R');
      val = t*trace(Xinv) - sum(log(x));
      S = V'*Xinv*V;   S2 = V'*Xinv^2*V;
      g = -t*diag(S2) - 1./x;
      H = t*2*(S.*S2) + diag(1./(x.^2));
      sol = [H ones(p,1); ones(1,p) 0] \ [-g;0];
      v = sol(1:p);
      fprime = v'*g;
      ntdecr = sqrt(max(-fprime,0));
      if (ntdecr < 1e-5), break; end;
      s = 1;
      newx = x+s*v;  newX = V*diag(newx)*V';
      [newR, flag] = chol(newX);
      while ((flag) | min(newx) < 0),
          s = BETA*s;
          newx = x+s*v;  newX = V*diag(newx)*V';
          [newR, flag] = chol(newX);
      end;
      newx = x+s*v;  newX = V*diag(newx)*V';
      [newR, flag] = chol(newX);
      newXinv = inv(newR)*inv(newR');
      newval = t*trace(newXinv) - sum(log(newx));
      while (newval > val + s*ALPHA*fprime)
          s = BETA*s;
          newx = x+s*v;  newX = V*diag(newx)*V';
          [newR, flag] = chol(newX);
          newXinv = inv(newR)*inv(newR');
          newval = t*trace(newXinv) - sum(log(newx));
      end;

      x = x + s*v;

   end;

   t = MU*t;
end;


figure(2)
% draw ellipsoid v'*W*v <= mu
W = inv(V*diag(x)*V')^2;
mu = diag(V'*W*V);
mu = mean(mu(ind));
angles = linspace(0,2*pi,noangles);
R = chol(W);  % W = R'*R
ellipsoid = sqrt(mu)*(R\[cos(angles); sin(angles)]);
d = plot(ellipsoid(1,:), ellipsoid(2,:), '--',0,0,'+');
set(d, 'Color', [0 0.5 0]);
set(d(2), 'MarkerFaceColor', [0 0.5 0]);
hold on

dot = plot(V(1,:),V(2,:),'o');
ind = find(x > 0.001);
dots = plot(V(1,ind),V(2,ind),'o');
set(dots,'MarkerFaceColor','blue');


xA = x;
disp('Nonzero lambda values for A design:');
for i=1:length(ind)
   text(V(1,ind(i)),V(2,ind(i)), ['l',int2str(ind(i))]);
   disp(['lambda(',int2str(ind(i)),') = ', num2str(x(ind(i)))]);
end;
%axis([-4.5 4.5 -4.5 4.5])
axis([-5 5 -5 5])
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
axis off
%print -deps Adesign.eps
hold off

keyboard



%%%%%%%%%% if (0) appears twice
% E-design

%
% minimize w
% s.t.     sum_i xi*vi*vi' >= w*I
%          x>= 0,  1'*x=1;
%

MU = 10;
ALPHA = 0.1;
BETA = 0.5;
TOL = 1e-5;

t = 1;
x = [(1/p)*ones(p,1); -1];   % w is last component of x

while ((p+2)/t>TOL)

   for i=1:100

      % minimize -t*w - logdet(X-w*I) - sum log xi
      % s.t.      1'*x = 1

      X = V*diag(x(1:p))*V' - x(p+1)*eye(2,2);
      R = chol(X);   % X = R'*R
      Xinv = inv(R)*inv(R');
      val = -t*x(p+1) - 2*sum(log(diag(R))) - sum(log(x(1:p)));
      S = V'*Xinv*V;   S2 = V'*(Xinv^2)*V;
      g = t*[zeros(p,1);-1] + [-diag(S); trace(Xinv)] - ...
          [1./x(1:p);  0] ;
      H =  [S.^2 -diag(S2); -diag(S2)' Xinv(:)'*Xinv(:)] ...
            + diag([1./(x(1:p).^2); 0]);
      sol = [H [ones(p,1);0] ; ones(1,p) 0 0] \ [-g;0];
      v = sol(1:(p+1));
      fprime = v'*g;
      ntdecr = sqrt(-fprime);
      if (ntdecr < 1e-5),
         W = Xinv/t;
         break;
      end;
      s = 1;
      newx = x+s*v;  newX = V*diag(newx(1:p))*V' - newx(p+1)*eye(2);
      [newR, flag] = chol(newX);
      while ((flag) | min(newx(1:p)) < 0),
          s = BETA*s;
          newx = x+s*v;
          newX = V*diag(newx(1:p))*V' - newx(p+1)*eye(2);
          [newR, flag] = chol(newX);
      end;
      newx = x+s*v;  newX = V*diag(newx(1:p))*V' - newx(p+1)*eye(2);
      [newR, flag] = chol(newX);
      newXinv = inv(newR)*inv(newR');
      newval = -t*newx(p+1) -  2*sum(log(diag(newR))) - ...
               sum(log(newx(1:p)));
      while (newval > val + s*ALPHA*fprime)
          s = BETA*s;
          newx = x+s*v;  newX = V*diag(newx(1:p))*V' - newx(p+1)*eye(2);
          [newR, flag] = chol(newX);
          newXinv = inv(newR)*inv(newR');
          newval = -t*newx(p+1) - 2*sum(log(diag(newR))) - ...
                   sum(log(newx(1:p)));
      end;

      x = x + s*v;

   end;

   t = MU*t;
end;


figure(3)
% draw ellipsoid v'*W*v <= mu
mu = diag(V'*W*V);
mu = mean(mu(ind));
angles = linspace(0,2*pi,noangles);
R = chol(W);  % W = R'*R
ellipsoid = sqrt(mu)*(R\[cos(angles); sin(angles)]);
d = plot(ellipsoid(1,:), ellipsoid(2,:), '--', 0, 0, '+');
set(d, 'Color', [0 0.5 0]);
set(d(2), 'MarkerFaceColor', [0 0.5 0]);
hold on

dot = plot(V(1,:),V(2,:),'o');
x = x(1:p);
ind = find(x > 0.001);
dots = plot(V(1,ind),V(2,ind),'o');
set(dots,'MarkerFaceColor','blue');


xE = x;
disp('Nonzero lambda values for E design:');
for i=1:length(ind)
   text(V(1,ind(i)),V(2,ind(i)), ['l',int2str(ind(i))]);
   disp(['lambda(',int2str(ind(i)),') = ', num2str(x(ind(i)))]);
end;
%axis([-4.5 4.5 -4.5 4.5])
axis([-5 5 -5 5])
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
axis off
%print -deps Edesign.eps
hold off

keyboard

%%%%%%%%%%%%%%%%%%


% E-design

%
% minimize w
% s.t.     sum_i xi*vi*vi' >= w*I
%          x>= 0,  1'*x=1;
%

MU = 10;
ALPHA = 0.1;
BETA = 0.5;
TOL = 1e-5;

t = 1;
x = [(1/p)*ones(p,1); -1];   % w is last component of x

while ((p+2)/t>TOL)

   for i=1:100

      % minimize -t*w - logdet(X-w*I) - sum log xi
      % s.t.      1'*x = 1

      X = V*diag(x(1:p))*V' - x(p+1)*eye(2,2);
      R = chol(X);   % X = R'*R
      Xinv = inv(R)*inv(R');
      val = -t*x(p+1) - 2*sum(log(diag(R))) - sum(log(x(1:p)));
      S = V'*Xinv*V;   S2 = V'*(Xinv^2)*V;
      g = t*[zeros(p,1);-1] + [-diag(S); trace(Xinv)] - ...
          [1./x(1:p);  0] ;
      H =  [S.^2 -diag(S2); -diag(S2)' Xinv(:)'*Xinv(:)] ...
            + diag([1./(x(1:p).^2); 0]);
      sol = [H [ones(p,1);0] ; ones(1,p) 0 0] \ [-g;0];
      v = sol(1:(p+1));
      fprime = v'*g;
      ntdecr = sqrt(-fprime);
      if (ntdecr < 1e-5),
         W = Xinv/t;
         break;
      end;
      s = 1;
      newx = x+s*v;  newX = V*diag(newx(1:p))*V' - newx(p+1)*eye(2);
      [newR, flag] = chol(newX);
      while ((flag) | min(newx(1:p)) < 0),
          s = BETA*s;
          newx = x+s*v;
          newX = V*diag(newx(1:p))*V' - newx(p+1)*eye(2);
          [newR, flag] = chol(newX);
      end;
      newx = x+s*v;  newX = V*diag(newx(1:p))*V' - newx(p+1)*eye(2);
      [newR, flag] = chol(newX);
      newXinv = inv(newR)*inv(newR');
      newval = -t*newx(p+1) -  2*sum(log(diag(newR))) - ...
               sum(log(newx(1:p)));
      while (newval > val + s*ALPHA*fprime)
          s = BETA*s;
          newx = x+s*v;  newX = V*diag(newx(1:p))*V' - newx(p+1)*eye(2);
          [newR, flag] = chol(newX);
          newXinv = inv(newR)*inv(newR');
          newval = -t*newx(p+1) - 2*sum(log(diag(newR))) - ...
                   sum(log(newx(1:p)));
      end;

      x = x + s*v;

   end;

   t = MU*t;
end;

% draw ellipsoid v'*W*v <= mu
mu = diag(V'*W*V);
mu = mean(mu(ind));
angles = linspace(0,2*pi,noangles);
R = chol(W);  % W = R'*R
ellipsoid = sqrt(mu)*(R\[cos(angles); sin(angles)]);
d = plot(ellipsoid(1,:), ellipsoid(2,:), '--', 0, 0, '+');
set(d, 'Color', [0 0.5 0]);
set(d(2), 'MarkerFaceColor', [0 0.5 0]);
axis([-5 5 -5 5])
hold on

dot = plot(V(1,:),V(2,:),'o');
x = x(1:p);
ind = find(x > 0.001);
dots = plot(V(1,ind),V(2,ind),'o');
set(dots,'MarkerFaceColor','blue');

xE = x;
disp('Nonzero lambda values for E design:');
for i=1:length(ind)
   text(V(1,ind(i)),V(2,ind(i)), ['l',int2str(ind(i))]);
   disp(['lambda(',int2str(ind(i)),') = ', num2str(x(ind(i)))]);
end;
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
axis off
% print -deps Edesign.eps
hold off


figure(4)

% confidence ellipsoids
% draw 90 percent confidence ellipsoid  for D design
W = V*diag(xD)*V';
angles = linspace(0,2*pi,noangles);
R = chol(W);  % W = R'*R
ellipsoid = sqrt(chi2inv(.9,3))*(R\[cos(angles); sin(angles)]);
plot(0,0,'ok',ellipsoid(1,:), ellipsoid(2,:), '-');
text(ellipsoid(1,1100),ellipsoid(2,1100),'D');
hold on

% draw 90 percent confidence ellipsoid  for A design
W = V*diag(xA)*V';
angles = linspace(0,2*pi,noangles);
R = chol(W);  % W = R'*R
ellipsoid = sqrt(chi2inv(.9,3))*(R\[cos(angles); sin(angles)]);
plot(0,0,'ok',ellipsoid(1,:), ellipsoid(2,:), '-');
text(ellipsoid(1,1),ellipsoid(2,1),'A');

% draw 90 percent confidence ellipsoid  for E design
W = V*diag(xE)*V';
angles = linspace(0,2*pi,noangles);
R = chol(W);  % W = R'*R
ellipsoid = sqrt(chi2inv(.9,3))*(R\[cos(angles); sin(angles)]);
d=plot(0,0,'ok',ellipsoid(1,:), ellipsoid(2,:), '-');
set(d,'Color',[0 0.5 0]);
text(ellipsoid(1,4000),ellipsoid(2,4000),'E');

% draw 90 percent confidence ellipsoid  for uniform design
W_u = inv(V*V'/p);
R = chol(W_u);  % W = R'*R
ellipsoid_u = sqrt(chi2inv(.9,3))*(R\[cos(angles); sin(angles)]);
plot(ellipsoid_u(1,:), ellipsoid_u(2,:), '--');
text(ellipsoid_u(1),ellipsoid_u(2),'U');
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
axis off
% print -deps confidence.eps
hold off

