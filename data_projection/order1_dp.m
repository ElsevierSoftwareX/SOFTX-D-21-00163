
%   ORDER 1 INITIAL CONDITION DATA PROJECTION
%
% We use the explicit 1st order expansion of the data projection algorethem that
% that is in Nicolsky et. al 2018
%
% Refernce = [Nicolsky et. al 2018]

function Phi = order1_dp(x)

  global eta0 u0

  eta0 = chebfun.interp1(x,eta0,'pchip');
  diff_eta0 = diff(eta0);
  u0 = chebfun.interp1(x,u0,'pchip');
  diff_u0 = diff(u0);

  ss = @(x) x + eta0(x);
  A = @(x) [0 1; ss(x) 0];
  B = [0 0; 1 0];

  D = @(x) (1+diff_eta0(x))*eye(2) + diff_u0(x)*A(x);

  phi0 = @(x) [u0(x) ; eta0(x)+(u0(x).^2)/2];
  phi0_prime = @(x) [diff_u0(x); diff_eta0(x)+u0(x)*diff_u0(x)];

  proj = @(x) phi0(x) + u0(x).*(diff_u0(x).*A(x)*inv(D(x))*B*phi0(x) -...
   B*phi0(x) - A(x)*inv(D(x))*phi0_prime(x));


  x_num = size(x);

  Phi = zeros(2, x_num(2));

  for i = 1:x_num(2)
    Phi(:, i) = proj(x(i));
  end

end
