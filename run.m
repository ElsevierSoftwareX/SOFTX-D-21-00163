
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

global eta_0 u_0 td l g           % physical variables
global t0 Tf x0 Xf t x            % physical variables
global x_res t_res                % grid resolution
global k la s n g

% plotting
WaveAnimation = 'on';
Eta_3D = 'off';
U_3D = 'off';
Runup = 'on';
L2Norm = 'on';
InitialConditions = 'off';
DataProjection = 'off';


% order of data projection
n = 3;

% doman of x and t
x0 = -2;
Xf = 10;
t0 = 0;
Tf = 3.5;

% setting resolution
x_res = 10000;      % number of points in the x domain - also used for s and k
t_res = 1000;       % number of points in the t domain - also used for lamda

% bathymetry and 'world' parameters
td = 1;           % slope of beach
g = 9.81;           % gravity acceleration
l = 1;              % arbitrary scaling parameter

% initializing arrays
[x,t,k,la,s] = Initialize();

% intial wave function parameters
H1 = 0.006;
H2 = 0.002;
c1 = 0.4;
c2 = 4.0;
x1 = 4.00;
x2 = 1.64;

eta_0 = H1*exp(-c1*(x-x1).^2)-H2*exp(-c2*(x-x2).^2);
%u_0 = 0.*x;
u_0 = -2*(sqrt(eta_0+td.*x)-sqrt(td.*x));

[eta_0, u_0] = IC(eta_0, u_0);

% checking intial conditions
STOP = CheckIC();

% computing solutions
[eta_analytic, u_analytic] = fast_hankel(n);
[eta_fvm, u_fvm] = run_num();

[ana, ana_u, num, num_u] = stat_norm(eta_analytic, eta_fvm, u_analytic, u_fvm);

plots = plotting(WaveAnimation, Eta_3D, U_3D, Runup, L2Norm, InitialConditions,...
  DataProjection);
