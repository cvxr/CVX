% Section 7.4.3: Probability bounds example with Voronoi diagram
% Boyd & Vandenberghe "Convex Optimization"
% Original version by Lieven Vandenberghe
% Adapted for CVX by Michael Grant, 2005/12/19
% Generates figures 7.5-7.8

randn('state',1);
rand('state',0);

% The vertices of the first voronoi set
Vs = [ 1 1 ; -1 2 ; -2 1 ; -2 -1 ; 0 -2 ; 1.5 -1 ; 1 1 ]';

% express first set as Ax <= b and compute the other constellation
% points
nopts = size( Vs, 2 );
m     = nopts - 1; 
As{1} = zeros(m,2);  
bs{1} = zeros(m,1);
Cs    = zeros( 2, nopts );
for i=1:m
   As{1}(i,:) = [-(Vs(2,i)-Vs(2,i+1))  Vs(1,i)-Vs(1,i+1)];
   bs{1}(i) = As{1}(i,:)*Vs(:,i);
   Cs(:,i+1)= 2* bs{1}(i)* As{1}(i,:)'/norm(As{1}(i,:))^2;
end;
figure(1)
plot([Vs(1,:) Vs(1,m)], [Vs(2,:), Vs(2,m)], '-', Cs(1,:), Cs(2,:), 'o');
hold on

% ellipses
for i=1 : nopts,
   noangles = 1000;
   text( Cs(1,i)+0.25, Cs(2,i)+0.25, [ 's_', int2str(i) ] );
   angles = linspace( 0, 2*pi, noangles );
%  R = sqrt(chi2inv(.5,2));  % 50% confidence ellipse
   R = 1;  % sigma = 1 ellipse
   ellipse = R * [ cos(angles) ; sin(angles) ] + Cs(:,i) * ones(1,noangles);
   plot( ellipse(1,:), ellipse(2,:), ':' );
end;

% Voronoi set around Cs(:,2)
As{2} = [-As{1}(1,:);
         (Cs(:,m+1) - Cs(:,2))';
         (Cs(:,3) - Cs(:,2))'];
bs{2} = [-bs{1}(1);
         As{2}(2,:)*(Cs(:,m+1)+Cs(:,2))/2;
         As{2}(3,:)*(Cs(:,3)+Cs(:,2))/2];
v = Vs(:,1);
dir = (Cs(:,m+1)+Cs(:,2))/2 - v;
plot( [v(1) v(1) + 5*dir(1)], [v(2)  v(2) + 5*dir(2)], '-');
v = Vs(:,2);
dir = (Cs(:,3)+Cs(:,2))/2 - v;
plot( [v(1) v(1) + 5*dir(1)], [v(2)  v(2) + 5*dir(2)], '-');

% Voronoi sets around Cs(:,i), i=3,4, .. nopts-1
for i=3:nopts-1
   As{i} = [-As{1}(i-1,:);
            (Cs(:,i-1) - Cs(:,i))';
            (Cs(:,i+1) - Cs(:,i))'];
   bs{i} = [-bs{1}(i-1);
            As{i}(2,:)*(Cs(:,i-1)+Cs(:,i))/2;
            As{i}(3,:)*(Cs(:,i+1)+Cs(:,i))/2];
   v = Vs(:,i-1);
   dir = (Cs(:,i-1)+Cs(:,i))/2 - v;
   plot( [v(1) v(1) + 4*dir(1)], [v(2)  v(2) + 4*dir(2)], '-');
   v = Vs(:,i);
   dir = (Cs(:,i+1)+Cs(:,i))/2 - v;
   plot( [v(1) v(1) + 4*dir(1)], [v(2)  v(2) + 4*dir(2)], '-');
end;

% Voronoi sets around Cs(:,nopts)
As{nopts} = [-As{1}(nopts-1,:);
             (Cs(:,nopts-1) - Cs(:,nopts))';
             (Cs(:,2) - Cs(:,nopts))'];
bs{nopts} = [-bs{1}(nopts-1);
             As{nopts}(2,:)*(Cs(:,nopts-1)+Cs(:,nopts))/2;
             As{nopts}(3,:)*(Cs(:,2)+Cs(:,nopts))/2];
v = Vs(:,nopts-1);
dir = (Cs(:,nopts-1)+Cs(:,nopts))/2 - v;
plot( [v(1) v(1) + 6*dir(1)], [v(2)  v(2) + 6*dir(2)], '-');
v = Vs(:,1);
dir = (Cs(:,2)+Cs(:,nopts))/2 - v;
plot( [v(1) v(1) + 6*dir(1)], [v(2)  v(2) + 6*dir(2)], '-');

axis('equal');
axis(5*[-1 1 -1 1]);
axis off
%print -deps chebbnds_example.eps

% Chebyshev lower bounds on probability of correct detection with
% sigma = 1 for set 1

% we calculate the chebyshev lower bound on
%  prob As{1}*(Cs(1) + v) <= bs{1}
A = As{1};
b = bs{1} - A*Cs(:,1);
[cd_cheb,P,q,r, X,lambda] = cheb(A,b,eye(2));

% to check the results, plot the ellipsoid
%  x'*P*x + 2*q'*x + r = 1

noangles = 500;
angles = linspace(0,2*pi, noangles);
pts = [cos(angles); sin(angles)];
ellipse = sqrt(1-r+q'*(P\q)) * P^(-1/2)*pts + ...
     (-P\q + Cs(:,1))*ones(1,noangles);
plot(ellipse(1,:), ellipse(2,:), 'r-');
dots= plot(X(1,:), X(2,:), 'ro');
set(dots,'MarkerFaceColor','red');
set(dots,'MarkerSize',4);

%print -deps chebbnds_example2.eps
hold off

%
% compute lower bounds vs sigma for a few sets
%

nosigmas = 100;
sigmas   = linspace( 0.001, 6.0, nosigmas )';
cd_cheb  = zeros( nosigmas, nopts );
cvxq     = cvx_quiet( true );
fprintf( 'Computing lower bounds' );
for k=1 : nosigmas,
    for i = 1 : 3,
        % Compute chebyshev lower bound for As{i}*(Cs(i)+v) <= bs{i}
        A = As{i};
        b = bs{i} - A*Cs(:,i);
        cd_cheb(k,i) = cheb( A, b, sigmas(k) * eye(2) );
   end;
   if rem( k, 10 ) == 0,
       fprintf( '.' );
   end
end;
fprintf( 'done.\n' );
cvx_quiet( cvxq );

figure(2)
plot(sqrt(sigmas(:,ones(1,3))), cd_cheb(:,[1 2 3]));
for i = 1 : 3,
   text( sqrt(sigmas(nosigmas/4)), cd_cheb(nosigmas/4,i), ['b',int2str(i)] );
end;
xlabel('x');
ylabel('y');
axis([0 2.5 0 1]);

%
% for central set, compute cheb lower bounds,  mc estimates,
% and chernoff bounds
%

nosigmas = 50;
sigmas   = linspace( 0.1, 0.5, nosigmas);
cd1      = zeros( 1, nosigmas );     % lower bounds for prob( x in C1 )
mc1      = zeros( 1, nosigmas );     % monte carlo estimates of prob( x in C1 )
cher1    = zeros( nopts-1, nosigmas ); % chernoff upper bounds on
                                     % Prob(x in Cj| s=s_1)
cvxq = cvx_quiet( true );
fprintf( 'Computing upper bounds, lower bounds, and Monte Carlo sims' );
for i = 1 : nosigmas,

   % calculate the chebyshev lower bound on
   % prob that As{1}*(Cs(1) + v) <= bs{1}

   A = As{1};  
   b = bs{1} - A*Cs(:,1);
   Sigma = sigmas(i) ^ 2 * eye(2);
   cd1(i) = cheb( A, b, Sigma );
   mc1(i) = montecarlo( A, b, Sigma, 10000 );

   % upper bounds on prob that Cs(1)+v in other sets
   % ie, As(j)*(Cs(1) + v ) <= bs(j)
   for j=2:nopts
       A = As{j};
       b = bs{j} - A*Cs(:,1);
       cher1( j-1, i ) = cher( A, b, Sigma );
   end;
   if rem( i, 5 ) == 0,
       fprintf( '.' );
   end
end;
cher1 = max( 1 - sum( cher1 ), 0 );
cvx_quiet( cvxq );
fprintf( 'done.\n' );

figure(4)
plot(sigmas, cher1, '-', sigmas, mc1, '--');
axis([0.2 0.5 0.9 1]);
xlabel('x'); 
ylabel('y');
%print -deps chernoff_example.eps
