% PDE solver in hodograph

% This solver uses a finite difference method to solve
% the NSWE in the (\sigma,\lambda) hodograph.


function [eta,u] = HodoSolve(Psi)

  global g td
  global t0 Tf t_res x_res numSig

  % setup parameters
  sig = linspace(0,1,numSig);           %sigma = 1 is boundary we care about
  lam = linspace(t0,Tf*sqrt(g),t_res);  %lambda, leave at 0-10

  % initialing
  dLam = lam(2)-lam(1);
  dSig = sig(2)-sig(1);
  courant = dLam/dSig;

  phi = zeros(t_res,numSig);
  psi = zeros(t_res,numSig);

  % boundary conditions
  phi(:, end) = Psi(1,:);
  psi(:, end) = Psi(2,:);
  %psi_bc = -0.003*sin(4*lam);

  % initial conditions
  %phi(1,:) = 0.00001*sin(sig);
  %psi(1,:) = H1*exp(-c1*(sig - x1).^2);



  fprintf('Please wait...\n');
  fprintf('Taking derivatives...\n')

% main loop and courant condition

A = zeros(2*numSig, 2*numSig);

r = dLam/(2*dSig);

%upper half of matrix   psi
for j = 2:numSig
  A(j,j) = 1;
  A(j,j+numSig-1) = -sig(j)*r;
  A(j,j+numSig) = dLam;
  A(j,j+numSig+1) = sig(j)*r;
end

%lower half of matrix phi
for j = numSig+1:2*numSig
  A(j, j-numSig-1) = -r;
  A(j, j-numSig+1) = r;
  A(j,j) = 1;
end

%implicit finite difference
for i=i:t_res

end

 %---------------PLOTTING------------------%

  figure(5)
  mesh(sig, lam, phi)
  title('Solution, $$\phi$$ in the hodograph plane','interpreter','latex')
  xlabel('$$\sigma$$','interpreter','latex')
  ylabel('$$\lambda$$','interpreter','latex')
  zlabel('$$\varphi$$','interpreter','latex')


  figure(6);
  mesh(sig, lam, psi);
  title('Solution, $$\psi$$ in the hodograph plane','interpreter','latex')
  xlabel('$$\sigma$$','interpreter','latex')
  ylabel('$$\lambda$$','interpreter','latex')
  zlabel('$$\psi$$','interpreter','latex')

  % figure(3);
  % mesh(courantList);

  fprintf('Graphing...\n')
  fprintf('Done.\n')

  %CG transform

  disp('    backwards CG transform and demensionalization... ');

  eta = zeros(t_res, numSig);
  u = zeros(t_res, numSig);
  xx = zeros(t_res, numSig);
  tt = zeros(t_res, numSig);

  for i=1:t_res
    %CG transform
    u(i,:) = phi(i,:);
    eta(i,:) = psi(i,:) - u(i,:).^2/2;
    xx(i,:) = sig - eta(i,:);
    tt(i,:) = u(i,:) + lam(i);

    %deminsionalizing
    u(i,:) = u(i,:)*sqrt(g*td);
    eta(i,:) = eta(i,:);
    tt(i,:) = tt(i,:)/sqrt(td*g);
  end

  %to display eta and u
    %figure(5);
    %mesh(tt,xx,eta);

    %figure(5);
    %mesh(tt,xx,u);

  disp('    scattered interpolation of eta and u... ');

  s_tt = reshape(tt, [t_res*numSig, 1]);
  s_xx = reshape(xx, [t_res*numSig, 1]);
  s_eta = reshape(eta, [t_res*numSig, 1]);
  s_u = reshape(u, [t_res*numSig, 1]);

  eta = scatteredInterpolant(s_xx, s_tt, s_eta);
  u = scatteredInterpolant(s_xx, s_tt, s_u);

end
