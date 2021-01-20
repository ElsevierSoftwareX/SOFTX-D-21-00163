function Phi = order_0_dp(x)
  global eta0 u0

  eta0 = chebfun.interp1(x,eta0,'pchip');
  u0 = chebfun.interp1(x,u0,'pchip');
  
  Phi = [u0(x) ; eta0(x)+(u0(x).^2)/2];

end
