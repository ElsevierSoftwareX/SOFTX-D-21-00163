%   ORDER N INITIAL CONDITION DATA PROJECTION


function Phi = order_n_dp(x)

    % Global variables:
    global eta_0 u_0
    global x_res t_res Xf g td n

    % Initial eta and u arrays at t=0
    eta0 = zeros(1,x_res);
    eta0(1:x_res) = eta_0(x);
    u0 = zeros(1,x_res);
    u0(1:x_res) = u_0(x);

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
    Phi0TCheb = u0gsCheb;                       % top component (phi)
    Phi0BCheb = eta0gsCheb+0.5*u0gsCheb.^2;     % bottom component (psi)
    dPhi0TdsCheb = diff(Phi0TCheb);             % derivatives with respect to sigma
    dPhi0BdsCheb = diff(Phi0BCheb);
    Phi0T = Phi0TCheb(gs0);                     % evaluating Psi0T at gamma(sigma)
    Phi0B = Phi0BCheb(gs0);                     % evaluating Psi0B at gamma(sigma)
    dPhi0Tds = dPhi0TdsCheb(gs0);               % evaluating dPsi0Tds at gamma(sigma)
    dPhi0Bds = dPhi0BdsCheb(gs0);               % evaluating dPsi0Bds at gamma(sigma)


    % Creating and inverting matrix D
    one = repelem(1,x_res);
    InvD = zeros(4,x_res);
    for i=1:x_res
      matrixD = [one(i),dPhi0Tds(i);s0(i)*dPhi0Tds(i),one(i)];
      InvD = inv(matrixD);
      InvLT(1,i) = InvD(1,1);                   % left-top element
      InvRT(1,i) = InvD(1,2);                   % right-top element
      InvLB(1,i) = InvD(2,1);                   % left-bottom element
      InvRB(1,i) = InvD(2,2);                   % right-bottom element
      InvD = [InvLT;InvRT;InvLB;InvRB];         % inverted matrix D (4xM matrix)
    end

    % empty arrays
    FkTStore = zeros(n,x_res);         % store F_k top
    FkBStore = zeros(n,x_res);         % store F_k bottom
    FkTds = zeros(1,x_res);            % derivative top
    FkBds = zeros(1,x_res);            % derivative bottom
    SumT = zeros(n,x_res);             % storing sum for each K value
    SumB = zeros(n,x_res);
    Phi_nT = zeros(1,x_res);           % solution arrays
    Phi_nB = zeros(1,x_res);

    % initializing FkT
    FkT = Phi0T;

    % main loop
    for i=1:n
      Term1 = ((Phi0T).^i)/factorial(i);        % non-recursive term
      FkTds = dPhi0Tds;                         % derivatives of F_k vector
      FkBds = dPhi0Bds;
      kphi = i*FkT;                             % k*phi
      Fk_NoInvD = [FkTds;kphi+FkBds];           % F_k vector before being mutiplied by InvD
      FkT = InvD(1,:).*Fk_NoInvD(1,:)+InvD(2,:).*Fk_NoInvD(2,:);        % mutiplying by InvD
      FkB = InvD(3,:).*Fk_NoInvD(1,:)+InvD(4,:).*Fk_NoInvD(2,:);
      FkTStore(i,:) = FkT;                      % storing FK values
      FkTBtore(i,:) = FkB;
      dPhi0TdsCheb = chebfun.interp1(s0(1,:),FkT(1,:),'pchip');   % refitting chebfuns to new F_k vector
      dPhi0BdsCheb = chebfun.interp1(s0(1,:),FkB(1,:),'pchip');
      dPhi0Tds = dPhi0TdsCheb(gs0);                               % evaluating F_k vector at gamma(sigma)
      dPhi0Bds = dPhi0BdsCheb(gs0);
      SumT(i,:) = Term1.*FkT;                   % storing for sumation
      SumB(i,:) = Term1.*FkB;
    end

    Phi_nT = Phi0T+sum(SumT);                   % computing nth order projection of phi
    Phi_nB = Phi0B+sum(SumB);                   % computing nth order projection of psi

    Phi = [Phi_nT; Phi_nB];

end
