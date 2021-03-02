function [x,t,k,la,s] = Initialize()

  global t0 Tf x0 Xf t x
  global t_res x_res g 

  % Initializing arrays
  x = linspace(x0,Xf,x_res);
  t = linspace(t0,Tf,t_res);
  k = linspace(0,70,x_res);
  la = linspace(t0,Tf*sqrt(g),t_res);
  s = linspace(0,Xf,x_res);

end
