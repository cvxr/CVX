% EX_8_5    One free point localization 
% Sec.8.7.3, Boyd & Vandenberghe "Convex Optimization"
% Joelle Skaf - 10/24/05
%
% K fixed points x_1,...,x_K in R^2 are given and the goal is to place 
% one additional point x such that the sum of the squares of the
% Euclidian distances to fixed points is minimized:
%           minimize    sum_{i=1}^K  ||x - x_i||^2 
% The optimal point is the average of the given fixed points

% Data generation 
n = 2;
K = 11;
randn('state',0);
P = randn(n,K);

% minimizing the sum of Euclidian distance 
fprintf(1,'Minimizing the sum of the squares the distances to fixed points...');

cvx_begin
    variable x(2) 
    minimize ( sum( square_pos( norms(x*ones(1,K) - P,2) ) ) )
cvx_end

fprintf(1,'Done! \n');

% Displaying results
disp('------------------------------------------------------------------');
disp('The optimal point location is: ');
disp(x);
disp('The average location of the fixed points is');
disp(sum(P,2)/K);
disp('They are the same as expected!');