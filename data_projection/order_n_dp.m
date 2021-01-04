%   ORDER N INITIAL CONDITION DATA PROJECTION
function Phi = order_n_dp(x)

    % Global variables:
    global eta_0 u_0
    global x_res t_res Xf g td n

    % Initial eta and u without dimensions
    [eta0, u0] = dimensionless(eta_0(x), u_0(x));

    % Computing sigma(x) at t=0
    s0 = x+eta0;

    % Interpolating inverse function of sigma -> gamma(sigma)
    gs0 = interp1(s0,x,x);
    gs0(1) = 0;

    % Interpolating u0(x) and eta0(x)
    u0gsCheb = chebfun.interp1(gs0, u0, 'pchip');        % evaluating at gamma(sigma)
    eta0gsCheb = chebfun.interp1(gs0, eta0, 'pchip');

    % Creating Phi0 vector
    Phi0phi = u0gsCheb;                                    %separating into top and bottom (phi and psi)
    Phi0psi = eta0gsCheb+0.5*u0gsCheb.^2;

    % taking individual derivatives with respect to sigma
    dPhi0phi_dsCheb = diff(u0gsCheb);
    dPhi0psi_dsCheb = diff(eta0gsCheb+0.5*u0gsCheb.^2);
    dPhi0phi_ds = dPhi0phi_dsCheb(s0);

    % Creating elements
    one = repelem(1,x_res);
    InvD = zeros(4,x_res);

    % Inverting the Matrix D
    for i=1:x_res
      matrixD = [one(i) dPhi0phi_ds(i);s0(i).*dPhi0phi_ds(i) one(i)];
      size(matrixD);
      InvD = inv(matrixD);
      InvLT(1,i) = InvD(1,1);                   % left-top element
      InvRT(1,i) = InvD(1,2);                   % right-top element
      InvLB(1,i) = InvD(2,1);                   % left-bottom element
      InvRB(1,i) = InvD(2,2);                   % right-bottom element
    end
    InvD = [InvLT;InvRT;InvLB;InvRB];         % assembling inverted matrix D (4xM matrix)

    % empty arrays
    SumT = zeros(n,x_res);            % storing sum for each K value
    SumB = zeros(n,x_res);            % storing sum for each K value
    Phi_n = zeros(2,x_res);           % solution arrays
    FkT = zeros(2,x_res);             % storing top of F_k vector for recursion
    FkB = zeros(2,x_res);             % storing bottom of F_k vector for recursion
    FkBds = zeros(1,x_res);           % derivative of bottom component of Fk vector
    FkT_NoInvD = zeros(1,x_res);
    FkB_NoInvD = zeros(1,x_res);
    % initializing FkT and FkB
    FkT(1,:) = Phi0phi(gs0);
    FkB(1,:) = Phi0psi(gs0);

    % main loop
    for i = 1:n
       phi_kfac = (((Phi0phi(gs0)).^i)/factorial(i));     % non-recursive term
       FkBCheb = FkB(1,:);
       FkBCheb = chebfun.interp1(s0,FkBCheb,'pchip');       % FkB to chebfun for derivative
       FkBdsCheb = -diff(FkBCheb);            % derivative of FkBCheb
       FkBds(1,:) = FkBCheb(s0);              % evaluating FkBCheb for discrete s0
       FkT_NoInvD = FkBds;                    % delta operator times FkT
       FkB_NoInvD = -2.*FkT;                  % delta operator times FkB
       FkT(2,:) = InvD(1,:).*FkT_NoInvD(1,:)+InvD(2,:).*FkT_NoInvD(1,:); % FkT for nth recursion
       FkB(2,:) = InvD(3,:).*FkB_NoInvD(1,:)+InvD(4,:).*FkB_NoInvD(1,:); % FkB for nth recursion
       SumT(i,:) = phi_kfac.*FkT(2,:);           % storing for later sumation
       SumB(i,:) = phi_kfac.*FkB(2,:);           % storing for later sumation
       FkT(2,:) = FkT(1,:);                      % preparing FkT array for next n
       FkB(2,:) = FkB(1,:);                      % preparing FkB array for next n
    end

    Phi_n(1,:) = Phi0phi(s0)+sum(SumT);
    Phi_n(2,:) = Phi0psi(s0)+sum(SumB);
    Phi = [Phi_n(1,:);Phi_n(2,:)];
  end
