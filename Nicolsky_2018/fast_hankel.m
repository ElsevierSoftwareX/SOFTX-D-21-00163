
%     FAST HANKEL SOLUTION TO 1 SPACIAL DEMENTION SWE
%
% We use an inverse hankel transfrom (not fast) for the computation of
% a(k) and b(k). Then we use a fast hankel transform for the compuatation of
% phi and psi. Finally, we use the CG transform to go back to original variables
% in (x,t) then return a scattered interpolant of eta and u.
%
% Refernce = [Nickolsky et. al. 2018]

function [eta u] = fast_hankel(n)

  disp('analytic solution... ');

  global x k la s                   %variables
  global x_res t_res Xf g td        %resolution
  global eta_0 u_0 eta0 u0
  global proj_phi proj_psi proj1_phi proj1_psi aofk bofk

  x_density = x_res/Xf;   % used for scaleing the inverse hankel transform

  % data_projection
  disp('    data_projection onto lambda = 0... ');
  proj = order_n_dp(s);
  proj1 = order1_dp(s);

  disp('    inverse hankel transform to compute a(k) and b(k)... ');

  % for the analytic solution to work properly both a and b need to be very
  % close to zero at the end of the domain of k
  a = 2*k.*ihat(proj(2, :) , sqrt(s), 2*k, 0)./x_density;
  b = -2*k.*ihat(proj(1, :) , sqrt(s), 2*k, 1)./x_density;

  disp('    fast hankel transform to compute psi and phi... ');

  for i=1:t_res  % for every lambda we do a fast hankel transform

    psi_freq = @(vk)  interp1(k,a,vk).*cos(la(i)*vk)+interp1(k,b,vk).*sin(la(i)*vk);
    phi_freq = @(vk)  interp1(k,a,vk).*sin(la(i)*vk)-interp1(k,b,vk).*cos(la(i)*vk);

    % fast hankel transform
    [psi(i, :) r_psi] = fht(psi_freq, 30, 7, 0, 20, 15);
    [phi(i,:) r_phi] = fht(phi_freq, 30, 7, 1, 20, 15);

    % scalling solution by 2*pi
    psi(i, :) = psi(i, :)./(2*pi);
    phi(i, :) = r_phi.^(1/2).*phi(i, :)./(2*pi);

  end

  r_size = size(psi); % same for phi and psi

  disp('    backwards CG transform and demensionalization... ');

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

  disp('    scattered interpolation of eta and u... ');

  s_tt = reshape(tt, [t_res*r_size(2), 1]);
  s_xx = reshape(xx, [t_res*r_size(2), 1]);
  s_eta = reshape(eta, [t_res*r_size(2), 1]);
  s_u = reshape(u, [t_res*r_size(2), 1]);

  eta = scatteredInterpolant(s_xx, s_tt, s_eta);
  u = scatteredInterpolant(s_xx, s_tt, s_u);

  % fitting chebfuns for plotting
  proj_phi = chebfun.interp1(s,proj(1,:),'pchip');
  proj_psi = chebfun.interp1(s,proj(2,:),'pchip');
  proj1_phi = chebfun.interp1(s,proj1(1,:),'pchip');
  proj1_psi = chebfun.interp1(s,proj1(2,:),'pchip');
  aofk = chebfun.interp1(k,a,'pchip');
  bofk = chebfun.interp1(k,b,'pchip');

end
