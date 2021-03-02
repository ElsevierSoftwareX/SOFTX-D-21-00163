
%Source Nicolsky 2018

%goal: add demensions to our independent and dependent variables

function [eta, u, x, t] = dimension(eta, u, x, t)

  global td g l

  eta = eta.*(l*td);
  u = u.*sqrt((g*td)/l);
  x = x.*l;
  t = t./(sqrt((g*td)/l));

end
