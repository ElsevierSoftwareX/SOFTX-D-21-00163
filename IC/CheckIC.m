function STOP = CheckIC()

  % global variables
  global eta_0 eta_prime u_0 u_prime td x_res x0 Xf

  % bathymetry Function (h)
  N  = x_res;							            % number of grid points
  x_check = linspace(x0, Xf, N+1)';			  % cell interfaces (the apostrophe is to transpose)
  dx = x_check(2) - x_check(1);                  	% spatial grid step
  xc = 0.5*(x_check(1:end-1) + x_check(2:end));  	% centers of cells
  h  = td.*xc;

  % evaluating functions for all x
  for i=1:length(xc)
    EtaCheck = eta_0(xc(i));
    EtaPrimeCheck = abs(eta_prime(xc(i)));
    uCheck = u_0(xc(i));

    if EtaCheck>=h
      STOP = 1;
      if STOP==1
        error('Error in Intial Condition: Reduce Amplitude of Eta.')
      end
    elseif EtaPrimeCheck>=1
      STOP = 2;
      if STOP==2
        error('Error in Intial Condition: Derivative of Eta is greater than 1 or -1.')
      end
    else
      STOP = 0;
      if STOP==0
        STOP=0;
      end
    end
  end
end
