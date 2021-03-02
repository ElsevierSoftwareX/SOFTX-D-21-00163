
%Source Nicolsky 2018

%goal: get rid of demensions such that the slope of a linear bathymetry is 1

function [eta, u, x, t] = dimensionless(eta, u, x, t)

  global td g l

  eta = eta./(l*td);
  u = u./(sqrt((g*td)/l));

  if nargin > 2
    x = x./l;
    t = t.*(sqrt((g*td)/l));
  end

end
