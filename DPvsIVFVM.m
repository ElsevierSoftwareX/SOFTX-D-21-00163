
clear
close all
format longE
%---------------RUNSCRIPT--------------------------%

addpath('IV_FVM/')
addpath('data_projection/')
addpath('FD_Hodograph/')

%%% Global variables:
global t0 Tf x0 Xf
global t_res x_res numSig
global eta_0 u_0 g td

% Bottom slope:
td = 10.0/10.0;

% initial time and final time (domain of t)
t0 = 0;
Tf = 0.55;

% domain of x
x0 = -0.1;
Xf = 10;

% resolution parameters
t_res = 200;
x_res = 2000;
numSig = 100;

%% Initial condition parameters:
H1 = 0.0006;
H2 = 0.0018;
c1 = 0.4444;
c2 = 4.0;
x1 = 4.1209+0.5;
x2 = 1.6384+0.5;

g = 9.81;

eta_0 = @(x) H1*exp(-c1*(x - x1).^2) - H2*exp(-c2*(x - x2).^2);
u_0   = @(x) 0;

[eta1, u1, eta_fvm, u_fvm] = BCpull();   % solution via FVM
Psi = order_n_BC_proj(eta1, u1, 4);   % Data Projection

  figure(3);
  plot(Psi(1,:));
  figure(4);
  plot(Psi(2,:));


[eta_hodo, u_hodo] = HodoSolve(Psi); % hodograph Solver

x_comp = round(x_res*(1-x0)/(Xf - x0), 0);

figure(7);
plot(eta_hodo(linspace(0,1,x_comp), linspace(0.699, 0.7, x_comp))), hold on;
plot(eta_fvm(1:x_comp, 199));


disp('post processing eta.....');

num = zeros(t_res, x_comp);
ana = zeros(t_res, x_comp);

t = linspace(t0, Tf, t_res);
x = linspace(x0, 1, x_comp);

for i = 1:t_res

  ana(i,:) = eta_hodo(linspace(x0, 1, x_comp), repmat(t(i), 1, x_comp));
  num(i,:) = eta_fvm(1:x_comp, i)';

  %making it zero on the other side of the shore
  max = 0;
  for j = 1:x_comp
    if num(i, j) + td*x(j) < 0
      num(i, j) = NaN;
      max = j;
    end
    if ana(i,j) + td*x(j) < 0
      ana(i,j) = NaN;
      max = j;
    end
  end
  stat_norm(i) = norm(ana(i, max+1:end) - num(i, max+1:end));
end

figure(8);
mesh(ana);

figure(9);
mesh(num);

figure(10);
plot(stat_norm);

disp('displaying .... ')

figure(11);
for i = 1:t_res
  disp(i);
  plot(x, num(i,:) ), hold on;
  plot(x, ana(i,:) ), hold on;
  plot(x, -td*x), hold off;
  axis([x0 1 -0.002 0.002]);
  pause(0.05);
end

figure(12);
subplot(4,1,1);
plot(x, num(1,:)), hold on;
plot(x, -td*x), hold on;
plot(x, ana(1,:)), hold off;
axis([x0 1 -0.002 0.002]);
xlabel('x');
ylabel('height');
title('t=0 (initial condition)');
legend('Numerical', 'Bathymetry', 'Analytical');

subplot(4,1,2);
plot(x, num(100,:)), hold on;
plot(x, -td*x), hold on;
plot(x, ana(100,:)), hold off;
axis([x0 1 -0.002 0.002]);
xlabel('x');
ylabel('height');
title('t=1 (run down)');
legend('Numerical', 'Bathymetry', 'Analytical');

subplot(4,1,3);
plot(x, num(200,:)), hold on;
plot(x, -td*x), hold on;
plot(x, ana(200,:)), hold off;
axis([x0 1 -0.002 0.002]);
xlabel('x');
ylabel('height');
title('t=2 (run down)');
legend('Numerical', 'Bathymetry', 'Analytical');

subplot(4,1,4);
plot(linspace(t0, Tf, t_res), stat_norm);
xlabel('t');
ylabel('L2 Norm')
title('L2 Norm at each time')
