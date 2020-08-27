
%   ORDER N BOUNDARY CONDITION DATA PROJECTION

function Psi = order_n_BC_proj(eta1, u1)

  global g t0 Tf t_res n

  % Parameters
  d=1;
  u1Dim = u1;                      % dimensional u
  eta1Dim = eta1;                  % dimensional eta
  tDim = linspace(t0,Tf,t_res);    % time array

  % Removing dimensions
  t = tDim/(sqrt(d/g));            % dimensionless time
  eta = eta1Dim/(d);               % dimensionless eta
  u = u1Dim/(sqrt(g*d));           % dimensionless u

  % Lambda
  Lambda = t-u;                     % Lambda array

  % Polynomial fitting to find gamma(Lambda)
  GammaCheb = chebfun.interp1(Lambda,t,'pchip');
  Gamma = GammaCheb(Lambda);

  % Computing eta at Gamma
  EtaGammaCheb = chebfun.interp1(Gamma,eta,'pchip');
  EtaGamma = EtaGammaCheb(Lambda);

  % Sigma
  SigmaCheb = EtaGammaCheb+1;

  % Creating Psi0 Vector
  Psi0T = Gamma-Lambda;                       % top component (phi)
  Psi0B = EtaGamma+0.5*(Gamma-Lambda).^2;     % bottom component (psi)

  % Interpolating Chebfuns for Psi0T and Psi0B
  Psi0TCheb = chebfun.interp1(Lambda,Psi0T,'pchip');
  Psi0BCheb = chebfun.interp1(Lambda,Psi0B,'pchip');

  % Differentiating Psi0T and Psi0B with respect to sigma
  dPsi0TdsCheb = diff(Psi0TCheb);
  dPsi0BdsCheb = diff(Psi0BCheb);
  dPsi0Tds = dPsi0TdsCheb(Lambda);    % evaluating at each time step
  dPsi0Bds = dPsi0BdsCheb(Lambda);

  % Compoinents of Matrix D
  dsdlCheb = diff(SigmaCheb);
  dsdl = dsdlCheb(Lambda);            % Derivative of sigma with respect to lambda
  sl = SigmaCheb(Lambda);             % Sigma(Lambda)
  MinusOne = repelem(-1,t_res);

  % creating and inverting matrix D
  InvD = zeros(4,t_res);
  for i=1:t_res
      MatrixD = [dsdl(i),MinusOne(i);-sl(i),dsdl(i)];
      InvD = inv(MatrixD);
      InvLT(1,i) = InvD(1,1);                   % left-top element
      InvRT(1,i) = InvD(1,2);                   % right-top element
      InvLB(1,i) = InvD(2,1);                   % left-bottom element
      InvRB(1,i) = InvD(2,2);                   % right-bottom element
      InvD = [InvLT;InvRT;InvLB;InvRB];         % inverted matrix D (4xM matrix

      % checking the deteriminant
      det0 = det(MatrixD);
      if det0 == 0
          disp('detD=0')
      end
  end

  % empty arrays
  FkTStore = zeros(n,t_res);          % store F_k top
  FkBStore = zeros(n,t_res);          % store F_k bottom
  FkTds = zeros(1,t_res);             % derivative top
  FkBds = zeros(1,t_res);             % derivative bottom
  SumT = zeros(n,t_res);              % storing sum for each K value
  SumB = zeros(n,t_res);
  Psi_nT = zeros(1,t_res);            % solution arrays
  Phi_nB = zeros(1,t_res);

  % initializing FkT
  FkT = Psi0T;

  % main loop
  for i=1:n
      Term1 = ((1-sl).^i)/factorial(i);     % non-recursive term
      FkTds = dPsi0Tds;                     % derivatives of F_k vector
      FkBds = dPsi0Bds;
      kphi = i*FkT;                         % k*phi
      Fk_NoInvD = [FkTds;kphi+FkBds];       % F_k vector before being mutiplied by InvD
      FkT = InvD(1,:).*Fk_NoInvD(1,:)+InvD(2,:).*Fk_NoInvD(2,:);   % mutiplying by InvD
      FkB = InvD(3,:).*Fk_NoInvD(1,:)+InvD(4,:).*Fk_NoInvD(2,:);
      FkTStore(i,:) = FkT;                  % storing FK values
      FkBStore(i,:) = FkB;
      dPsi0TdsCheb = chebfun.interp1(Lambda(1,:),FkT(1,:),'pchip');  % refitting chebfuns to new F_k vector
      dPsi0BdsCheb = chebfun.interp1(Lambda(1,:),FkB(1,:),'pchip');
      dPhi0Tds = dPsi0TdsCheb(Lambda);                               % evaluating F_k vector
      dPhi0Bds = dPsi0BdsCheb(Lambda);
      SumT(i,:) = Term1.*FkT;               % storing for sumation
      SumB(i,:) = Term1.*FkB;
  end

  Psi_nT = Psi0T+sum(SumT);         % computing nth order projection of phi
  Psi_nB = Psi0B+sum(SumB);         % computing nth order projection of psi

  Psi = [Psi_nT; Psi_nB];

end
