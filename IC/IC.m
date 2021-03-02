function [eta_0, u_0] = IC(eta_0, u_0);

  global x x_res eta_prime u_prime

  eta_0 = chebfun.interp1(x,eta_0,'pchip');
  eta_prime = diff(eta_0);
  for i=1:x_res
    u_0(imag(u_0)~=0) = 0;
  end

  u_0 = chebfun.interp1(x,u_0,'pchip');
  u_prime = diff(u_0);
end
