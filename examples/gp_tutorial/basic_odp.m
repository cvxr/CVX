% Optimal doping profile optimization
% Boyd, Kim, Vandenberghe, and Hassibi, "A tutorial on geometric programming"
% Joshi, Boyd, and Dutton, "Optimal doping profiles via geometric programming"
% Written for CVX by Almir Mutapcic 02/08/06
% (a figure is generated)
%
% Determines the optimal doping profile that minimizes base transit
% time in a (homojunction) bipolar junction transistor.
% This problem can be posed as a GP:
%
%   minimize   tau_B
%       s.t.   Nmin <= v <= Nmax
%              y_(i+1) + v_i^const1 <= y_i
%              w_(i+1) + v_i^const2 <= w_i, etc...
%
% where variables are v_i, y_i, and w_i.

% discretization size
M = 50;
% The old version of this example took a long time to build with M=1000,
% because it used for loops to build the y/v/w constraints. This version
% vectorizes these loops, so 1000 points are handled in a couple of 
% seconds. I do recommend SeDuMi or MOSEK for this, however.
% M = 1000; 

% problem constants
g1 = 0.42;
g2 = 0.69;
Nmax = 5*10^18;
Nmin = 5*10^16;
Nref = 10^17;
Dn0 = 20.72;
ni0 = 1.4*(10^10);
WB = 10^(-5);
C =  WB^2/((M^2)*(Nref^g1)*Dn0);

% exponent powers
pwi = g2 -1;
pwj = 1+g1-g2;

cvx_begin gp
    variables v(M) y(M) w(M)
    % objective function is the base transmit time
    tau_B = C*w(1);
    minimize( tau_B )
    subject to
        Nmin <= v <= Nmax; %#ok
        y >= [ y(2:end) ; 0 ] + v .^ pwj; %#ok
        w >= [ w(2:end) ; 0 ] + y .* v .^ pwi; %#ok
cvx_end

% plot the basic optimal doping profile
figure, clf
nbw = 0:1/M:1-1/M;
semilogy(nbw,v,'LineWidth',2);
axis([0 1 1e16 1e19]);
xlabel('base');
ylabel('doping');
text(0,Nmin,'Nmin ', 'HorizontalAlignment','right');
text(0,Nmax,'Nmax ', 'HorizontalAlignment','right');
disp('Optimal doping profile is plotted.')
