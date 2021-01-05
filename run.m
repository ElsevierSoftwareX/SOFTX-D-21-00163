
%   GENERAL INITIAL CONDITION ANALYTIC VS NUMERICAL SOLUTION COMPARISON
%
% Works with a given IC on a given domain. Produces a comparison of the analytic
% solution
%
% Notes:
%    a(k) and b(k) in "Nicolsky_2018/fast_hankel" must converge to zero
%       before the end of the k domain. Please adjust both the resolution
%       and domain for this to effectively work
%

clear
close all
format longE

global eta_0 eta_prime u_0 u_prime td l g t0 Tf x0 Xf t % physical variables
global x_res t_res                                      % grid resolution
global x k la s n                                       % variables
global Initial_Eta Initial_u Data_Projection a_of_k     %
global b_of_k Time_Frames L2_Norm Wave_Animation
global num ana stat_norm stat_norm_max

% plotting
Initial_Eta = 'off';
Initial_u = 'off';
Data_Projection = 'off';
a_of_k = 'off';
b_of_k = 'off';
L2_Norm = 'on';
Wave_Animation = 'off';


% order of data projection
n = 5;

% doman of x and t
x0 = -2;
Xf = 10;
t0 = 0;
Tf = 5;

% setting resolution
x_res = 500;     % number of points in the x domain - also used for s and k
t_res = 201;      % number of points in the t domain - also used for lamda

% intial wave function parameters
H1 = 0.006;
H2 = 0.018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209;
x2 = 1.6384;

% defining initial conditions (as chebfuns)
eta_0 = chebfun(@(x) H1*exp(-c1*(x-x1).^2)-H2*exp(-c2*(x-x2).^2), [x0, Xf]);
u_0   = chebfun(@(x) -0.05*sin(3*x)*exp(-0.5*(x-5)^2), [x0 Xf]);
eta_prime = diff(eta_0);
u_prime = diff(u_0);


% bathymetry and 'wolrd' parameters
td = 10.0/10.0;     % slope of beach
g = 9.81;           % gravity acceleration
l = 1;              % abritrary scalling parameter

% initializing arrays
x = linspace(x0,Xf,x_res);
t = linspace(t0,Tf,t_res);
k = linspace(0,70,x_res);
la = linspace(t0,Tf*sqrt(g),t_res);
s = linspace(0,Xf,x_res);

% checking intial conditions
STOP = CheckIC();

% computing solutions
[eta_analytic u_analytic] = fast_hankel(n);
[eta_fvm u_fvm] = run_num();

%post processing
num = zeros(t_res, x_res);
ana = zeros(t_res, x_res);

for i = 1:t_res

  ana(i,:) = eta_analytic(x, repmat(t(i), 1, x_res));
  num(i,:) = eta_fvm(:, i)';

  % keeping wave above beach in plotting
  maximum = 0;
  for j = 1:x_res
    if num(i, j) + td*x(j) < 0
      num(i, j) = NaN;
      maximum = j;
    end
    if ana(i,j) + td*x(j) < 0
      ana(i,j) = NaN;
      maximum = j;
    end
  end
  stat_norm(i) = norm(ana(i, maximum+1:end) - num(i, maximum+1:end));
end

stat_norm_max = find(stat_norm == max(stat_norm));

plots = plotting();
