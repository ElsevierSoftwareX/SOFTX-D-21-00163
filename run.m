
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

clear all
close all
format longE

global eta_0 u_0 td l g               % physical variables
global t0 Tf x0 Xf t                                    % physical variables
global x_res t_res                                      % grid resolution
global x k la s n                                       % variables
global eta_analytic eta_fvm max_ana_x                          % post processing/stats

% plotting
Initial_Eta = 'off';
Initial_u = 'off';
Data_Projection = 'on';
a_of_k = 'off';
b_of_k = 'off';
L2_Norm = 'off';
Wave_Animation = 'off';


% order of data projection
n = 1;

% doman of x and t
x0 = -2;
Xf = 10;
t0 = 0;
Tf = 3;

% setting resolution
x_res = 10000;      % number of points in the x domain - also used for s and k
t_res = 10000;       % number of points in the t domain - also used for lamda

% bathymetry and 'world' parameters
td = 1.0;           % slope of beach
g = 9.81;           % gravity acceleration
l = 1;              % arbitrary scaling parameter

% initializing arrays
x = linspace(x0,Xf,x_res);
t = linspace(t0,Tf,t_res);
k = linspace(0,70,x_res);
la = linspace(t0,Tf*sqrt(g),t_res);
s = linspace(0,Xf,x_res);


% intial wave function parameters
H1 = 0.006;
H2 = 0.018;
c1 = 0.44;
c2 = 4.0;
x1 = 4.12;
x2 = 1.64;
eta_0 = H1*exp(-c1*(x-x1).^2)-H2*exp(-c2*(x-x2).^2);
u_0 = 0.*x;
%u_0 = -2*(sqrt(eta_0+td.*x)-sqrt(td.*x));

[eta_0, u_0] = IC(eta_0, u_0);

% checking intial conditions
STOP = CheckIC();

% computing solutions
[eta_analytic, u_analytic] = fast_hankel(n);
[eta_fvm, u_fvm] = run_num();

stats = stat_norm();

plots = plotting(Initial_Eta, Initial_u, Data_Projection,...
  a_of_k, b_of_k, L2_Norm, Wave_Animation);
