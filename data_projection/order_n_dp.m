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
    gs0Cheb = chebfun.interp1(s0(1,:),x(1,:), 'pchip');
    gs0 = gs0Cheb(x(1,:));      % evaluating at x

    % Interpolating u0(x) and eta0(x)
    u0Cheb = chebfun.interp1(x(1,:),u0(1,:),'pchip');
    eta0Cheb = chebfun.interp1(x(1,:),eta0(1,:),'pchip');
    u0gs = u0Cheb(gs0);        % evaluating at gamma(sigma)
    eta0gs = eta0Cheb(gs0);

    % Interpolating sigma(u0gs) and sigma(eta0gs)
    eta0gsCheb = chebfun.interp1(s0(1,:),eta0gs(1,:),'pchip');
    u0gsCheb = chebfun.interp1(s0(1,:),u0gs(1,:),'pchip');

    % Creating Phi0 vector
    Phi0Cheb = [u0gsCheb eta0gsCheb+0.5*u0gsCheb.^2];     % Phi0 with chebfuns
    dPhi0dsCheb = diff(Phi0Cheb);                         % derivatives with respect to sigma
    Phi0 = [Phi0Cheb(gs0, 1); Phi0Cheb(gs0, 2)];          % evaluating Psi0 at gamma(sigma)
    dPhi0ds = [dPhi0dsCheb(gs0, 1); dPhi0dsCheb(gs0, 2)]; % evaluating dPsi0Tds at gamma(sigma)

    % Creating and inverting matrix D
    one = repelem(1,x_res);
    InvD = zeros(4,x_res);
    for i=1:x_res
      matrixD = [one(i),dPhi0ds(1,i);s0(i)*dPhi0ds(2,i),one(i)];
      InvD = inv(matrixD);
      InvLT(1,i) = InvD(1,1);                   % left-top element
      InvRT(1,i) = InvD(1,2);                   % right-top element
      InvLB(1,i) = InvD(2,1);                   % left-bottom element
      InvRB(1,i) = InvD(2,2);                   % right-bottom element
      InvD = [InvLT;InvRT;InvLB;InvRB];         % inverted matrix D (4xM matrix)
    end

    % empty arrays
    Fkds = zeros(1,x_res);            % derivative
    SumT = zeros(n,x_res);            % storing sum for each K value
    SumB = zeros(n,x_res);            % storing sum for each K value
    Phi_n = zeros(1,x_res);           % solution arrays

    % initializing FkT
    FkT = Phi0(1,:);

    % main loop
    for i=1:n
      Term1 = ((Phi0(1,:)).^i)/factorial(i);        % non-recursive term
      Fkds = dPhi0ds;                               % derivatives of F_k vector
      kphi = i*FkT;                                 % k*phi
      Fk_NoInvD = [Fkds(1,:);kphi+Fkds(2,:)];       % F_k vector before being mutiplied by InvD
      Fk = [InvD(1,:).*Fk_NoInvD(1,:)+InvD(2,:).*Fk_NoInvD(2,:); InvD(3,:).*Fk_NoInvD(1,:)+InvD(4,:).*Fk_NoInvD(2,:)];       % mutiplying by InvD
      dPhi0TdsCheb = chebfun.interp1(s0(1,:),Fk(1,:),'pchip');
      dPhi0BdsCheb = chebfun.interp1(s0(1,:),Fk(2,:),'pchip');  % refitting chebfuns to new F_k vector
      dPhi0ds = [dPhi0TdsCheb(gs0); dPhi0BdsCheb(gs0)];                               % evaluating F_k vector at gamma(sigma)
      SumT(i,:) = Term1.*Fk(1,:);                   % storing for sumation
      SumB(i,:) = Term1.*Fk(2,:);                   % storing for sumation
    end

    Phi = [Phi0(1,:)+sum(SumT); Phi0(2,:)+sum(SumB)];                   % computing nth order projection of Phi

end
