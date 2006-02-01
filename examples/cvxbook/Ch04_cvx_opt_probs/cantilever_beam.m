% Exercise 4.31: Design of a cantilever beam (non-recursive convex GP formulation)
% (For a detailed explanation see section 4.5.4, pp. 163-165)
% Boyd & Vandenberghe "Convex Optimization"
% Almir Mutapcic - 01/30/06
% (a figure is generated)
%
% We have a segmented cantilever beam with N segments. Each segment
% has a unit length and variable width and height (rectangular profile).
% The goal is minimize the total volume of the beam, over all segment
% widths w_i and heights h_i, subject to constraints on aspect ratios,
% maximum allowable stress in the material, vertical deflection y, etc.
%
% The problem can be posed as a geometric program (posynomial form)
%     minimize    sum( w_i* h_i)
%         s.t.    w_min <= w_i <= w_max,       for all i = 1,...,N
%                 h_min <= h_i <= h_max
%                 S_min <= h_i/w_i <= S_max
%                 6*i*F/(w_i*h_i^2) <= sigma_max
%                 6*F/(E*w_i*h_i^3) == d_i
%                 (2*i - 1)*d_i + v_(i+1) <= v_i
%                 (i - 1/3)*d_i + v_(i+1) + y_(i+1) <= y_i
%                 y_1 <= y_max
%
% with variables w_i, h_i, d_i, (i = 1,...,N) and v_i, y_i (i = 1,...,N+1).
% (Consult the book for other definitions and a recursive formulation of
% this problem.)

% data generation 
N = 4;

% constants
wmin = .1; wmax = 100;
hmin = .1; hmax = 6;
Smin = 1/5; Smax = 5;
sigma_max = 1;
ymax = 10;
E = 1; F = 1;

% formulating the problem as a GP in the convex form
cvx_begin
  variables wlog(N) hlog(N) vlog(N+1) ylog(N+1) dlog(N)

  minimize ( logsumexp_sdp(wlog + hlog) )
  subject to

  % constraints on the rectangular profile variables
  wlog >= log(wmin)
  wlog <= log(wmax)
  hlog >= log(hmin)
  hlog <= log(hmax)
  hlog - wlog >= log(Smin)
  hlog - wlog <= log(Smax)

  % maximum stress constraint
  log(6*F*[1:N]') - (wlog + 2*hlog) <= log(sigma_max)

  % force and deflection constraints
  log(6*F) - (log(E) + wlog + 3*hlog) == dlog
  for i = 1:N
    logsumexp_sdp( [log(2*i-1)+dlog(i) vlog(i+1)] ) <= vlog(i)
    logsumexp_sdp( [log(i-1/3)+dlog(i) vlog(i+1) ylog(i+1)] ) <= ylog(i)
  end
  ylog(1) <= log(ymax);
cvx_end

w = exp(wlog);
h = exp(hlog);

% display results 
disp('The optimal widths and heights are: ');
w, h
fprintf(1,'The optimal minimum volume of the beam is %3.4f\n', sum(w.*h))

% plot the 3D model of the optimal cantilever beam
close all;
cantilever_beam_plot([h; w])
