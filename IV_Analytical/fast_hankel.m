
%     FAST HANKEL SOLUTION TO 1 SPATIAL DIMENSION SWE
%
% We use an inverse hankel transfrom (not fast) for the computation of
% a(k) and b(k). Then we use a fast hankel transform for the compuatation of
% phi and psi. Finally, we use the CG transform to go back to original variables
% in (x,t) then return a scattered interpolant of eta and u.
%
% Refernce = [Nickolsky et. al. 2018]

function [eta, u] = fast_hankel(n)

  tStart = tic;
  fprintf('\n')
  disp('Analytic Solution:');

  % global variables
  global Xf x_res t_res k la s proj0 proj1 projn

  x_density = x_res/Xf;   % used for scaling the inverse hankel transform

  % data_projection
  display_hankel = fprintf('   Projecting IC onto lambda = 0.');
  projn = order_n_dp(s);
  proj0 = order_0_dp(s);
  proj1 = order1_dp(s);

  fprintf(repmat('\b',1,display_hankel));
  display_hankel = fprintf('   Inverse hankel transfrom to compute a(k) and b(k).');

  % for the analytic solution to work properly both a and b need to be very
  % close to zero at the end of the domain of k

  if n == 0
    a = 2*k.*ihat(proj0(2, :) , sqrt(s), 2*k, 0)./x_density;
    b = -2*k.*ihat(proj0(1, :) , sqrt(s), 2*k, 1)./x_density;
  elseif n == 1
    a = 2*k.*ihat(proj1(2, :) , sqrt(s), 2*k, 0)./x_density;
    b = -2*k.*ihat(proj1(1, :) , sqrt(s), 2*k, 1)./x_density;
  else
    a = 2*k.*ihat(projn(2, :) , sqrt(s), 2*k, 0)./x_density;
    b = -2*k.*ihat(projn(1, :) , sqrt(s), 2*k, 1)./x_density;
  end


  fprintf(repmat('\b',1,display_hankel));
  display_hankel = fprintf('   Fast hankel transform to compute psi and phi.');

  for i=1:t_res  % for every lambda we do a fast hankel transform

    psi_freq = @(vk)  interp1(k,a,vk).*cos(la(i)*vk)+interp1(k,b,vk).*sin(la(i)*vk);
    phi_freq = @(vk)  interp1(k,a,vk).*sin(la(i)*vk)-interp1(k,b,vk).*cos(la(i)*vk);

    % fast hankel transform
    [psi(i, :), r_psi] = fht(psi_freq, 30, 7, 0, 20, 15);
    [phi(i,:), r_phi] = fht(phi_freq, 30, 7, 1, 20, 15);

    % scaling solution by 2*pi
    psi(i, :) = psi(i, :)./(2*pi);
    phi(i, :) = r_phi.^(1/2).*phi(i, :)./(2*pi);

  end

  r_size = size(psi); % same for phi and psi

  fprintf(repmat('\b',1,display_hankel));
  display_hankel = fprintf('   Backwards CG transform and dimensionalization.');

  eta = zeros(t_res, r_size(2));
  u = zeros(t_res, r_size(2));
  xx = zeros(t_res, r_size(2));
  tt = zeros(t_res, r_size(2));

  for i=1:r_size(2)
    %CG transform
    u(:,i) = phi(:,i);
    eta(:,i) = psi(:,i) - u(:,i).^2/2;
    xx(:,i) = r_phi(i).^2./4 - eta(:,i);
    tt(:,i) = u(:,i) + la';
  end

  %deminsionalizing
  [eta, u, xx, tt] = dimension(eta, u, xx, tt);

  fprintf(repmat('\b',1,display_hankel));
  display_hankel = fprintf('   Scattered Interpolation of eta and u.');

  s_tt = reshape(tt, [t_res*r_size(2), 1]);
  s_xx = reshape(xx, [t_res*r_size(2), 1]);
  s_eta = reshape(eta, [t_res*r_size(2), 1]);
  s_u = reshape(u, [t_res*r_size(2), 1]);

  eta = scatteredInterpolant(s_xx, s_tt, s_eta);
  u = scatteredInterpolant(s_xx, s_tt, s_u);

  tEnd = toc(tStart);
  fprintf(repmat('\b',1,display_hankel));
  display_hankel = fprintf('   Computation time: %0.4f \n \n', tEnd);
end
