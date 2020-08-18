function STOP = CheckIC()

  % Global variables
  global eta_0 eta_prime u_0 u_prime td x_res x0 Xf

  % Bathymetry Function (h)
  N  = x_res;							            % number of grid points
  x  = linspace(x0, Xf, N+1)';			  % cell interfaces (the apostrophe is to transpose)
  dx = x(2) - x(1);                  	% spatial grid step
  xc = 0.5*(x(1:end-1) + x(2:end));  	% centers of cells
  h  = td.*xc;

  % evaluating functions for all x
  for i=1:length(xc)
    EtaCheck = eta_0(xc(i));
    EtaPrimeCheck = abs(eta_prime(xc(i)));
    uCheck = u_0(xc(i));

    if EtaCheck>=h
      STOP = 1;
    elseif EtaPrimeCheck>=1
      STOP = 2;
    elseif uCheck>=0.05
      STOP = 3;
    else
      STOP = 0;
    end

    if STOP==1
      error('Error in Intial Condition: Reduce Amplitude of Eta.')
    end
    if STOP==2
      error('Error in Intial Condition: Derivative of Eta is greater than 1 or -1.')
    end
    if STOP==3
      error('Error in Intial Condition: Reduce initial velocity.')
    end
    if STOP==0
      STOP=0;
    end
  end
end
