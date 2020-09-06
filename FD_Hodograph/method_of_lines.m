%---PDE solver in Hodograph---------%

% This solver uses a finite difference method to solve
% the NSWE in the (\sigma,\lambda) hodograph.


function [eta,u] = HodoSolve(Psi)

  global g td
  global t0 Tf t_res x_res numSig

  %--------------SETUP PARAMTERS------------------%

  sig = linspace(0,1,numSig); %sigma = 1 is boundary we care about
  lam = linspace(t0,Tf*sqrt(g),t_res); %lambda, leave at 0-10

  %-----------INITIALIZING---------%

  dLam = lam(2)-lam(1);
  dSig = sig(2)-sig(1);
  courant = dLam/dSig;

  phi = zeros(t_res,numSig);
  psi = zeros(t_res,numSig);

  %-------------BOUNDARY CONDITIONS------------%

  phi(:, end) = Psi(1,:);
  psi(:, end) = Psi(2,:);
  %psi_bc = -0.003*sin(4*lam);

  %------------INITIAL CONDITIONS---------------%

  %phi(1,:) = 0.00001*sin(sig);
  %psi(1,:) = H1*exp(-c1*(sig - x1).^2);



  fprintf('Please wait...\n');
  fprintf('Taking derivatives...\n')

%-----------MAIN LOOP AND COURANT CONDITION------------%

options = odeset('AbsTol', 1e-4, 'RelTol', 1e-4, 'OutputFcn', @odetpbar, 'MaxStep', 1.0);
 = ode23(@pde1, [t0 Tf], w0, options);


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


  %define the conservatiuon of mass with description of space but not in time
  function dpsi = pde1(lambda, psi, phi)

    %for a given lambda this is the psi and phi

    for i=1:(numSig-1)
      dpsi(i) = -(phi(i) - phi(i+1))./dSig - phi(i);
    end

  end

  %difine the conservatiuon of momentum with descriatipon of space but not in time
  function dphi = pde2(lambda, psi, phi)

    for i=1:(numSig-1)
      dphi(i) = -(psi(i) - psi(i+1))./dSig;
    end

  end

end
