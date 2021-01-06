function plots = plotting()

  % function that produces all plots
  global Initial_Eta Initial_u Data_Projection a_of_k...
    b_of_k L2_Norm Wave_Animation eta_0 u_0 s n proj_phi...
    proj_psi proj1_phi proj1_psi k aofk bofk x td x0 Xf t0...
    Tf t t_res stat_norm num ana stat_norm_max eta0 u0

  cnt = 1;   % figure counter

  % initial eta plot
  if contains(Initial_Eta,"on")
    figure(cnt)
    plot(eta_0)
    title('Initial Wave Height ($\eta_0$)','Interpreter','latex')
    xlabel('$x$','Interpreter','latex')
    ylabel('$\eta_0$','Interpreter','latex')
    cnt = cnt+1;
  else
    ;
  end

  % initial velocity plot
  if contains(Initial_u,"on")
    figure(cnt)
    plot(u_0)
    title('Initial Wave Velocity ($u_0$)','Interpreter','latex')
    xlabel('$x$','Interpreter','latex')
    ylabel('$u_0$','Interpreter','latex')
    cnt = cnt+1;
  else
    ;
  end

  % plot of data projection
  if contains(Data_Projection,"on")
    figure(cnt);
    plot(proj_phi), hold on;
    plot(proj1_phi), hold on;
    plot(x, u0(x));
    title('Projection of $\varphi_0$', 'Interpreter','latex');
    xlabel('$\sigma$', 'Interpreter','latex');
    ylabel('$\varphi$', 'Interpreter','latex');
    legend({sprintf('%dth order projection',n)','$1^{st}$ order projection',...
      '$\varphi_0$'},'Location','southwest','Interpreter','latex');
    cnt = cnt+1;

    figure(cnt)
    plot(proj_psi), hold on;
    plot(proj1_psi), hold on;
    plot(x, eta0(x)+u0(x).^2/2), hold on;
    title('Projection of $\psi_0$', 'Interpreter','latex');
    xlabel('$\sigma$', 'Interpreter','latex');
    ylabel('$\psi_n$', 'Interpreter','latex');
    legend({sprintf('%dth order projection', n),'$1^{st}$ order projection',...
     '$\psi_0$'},'Location','southwest','Interpreter','latex');
     cnt = cnt+1;
   else
     ;
   end

   % plot of a(k)
   if contains(a_of_k,"on")
     figure(cnt)
     plot(aofk)
     title('$a(k)$', 'Interpreter','latex');
     xlabel('$k$', 'Interpreter','latex');
     ylabel('$a$', 'Interpreter','latex');
     cnt = cnt+1;
   else
     ;
   end

   % plot of b(k)
   if contains(b_of_k,"on")
     figure(cnt)
     plot(bofk)
     title('$b(k)$', 'Interpreter','latex');
     xlabel('$k$', 'Interpreter','latex');
     ylabel('$b$', 'Interpreter','latex');
     cnt = cnt+1;
   else
     ;
   end

   % L2 norm plot
   if contains(L2_Norm,"on")
     Th = Tf/2;
     TL2Max = Tf/t_res*stat_norm_max;   % time with max l2 norm

     figure(cnt)
     subplot(4,1,1);
     plot(x,num(1,:)), hold on;
     plot(x,ana(1,:)), hold on;
     plot(x,-td*x), hold off;
     axis([-0.2 6 -0.03 0.04]);
     title('Initial Condition ($t=0$)', 'Interpreter','latex');
     xlabel('$x$', 'Interpreter','latex');
     ylabel('$\eta$', 'Interpreter','latex');
     legend('Numerical', 'Analytical','Bathymetry');

     subplot(4,1,2);
     plot(x,num(stat_norm_max,:)), hold on;
     plot(x,ana(stat_norm_max,:)), hold on;
     plot(x,-td*x), hold off;
     axis([-0.2 6 -0.03 0.04]);
     title(sprintf('Time with maximum L2 norm ($t=%0.2f$)',TL2Max), 'Interpreter','latex');
     xlabel('$x$', 'Interpreter','latex');
     ylabel('$\eta$', 'Interpreter','latex');
     legend('Numerical', 'Analytical','Bathymetry');

     subplot(4,1,3);
     plot(x,num(floor(0.5*t_res),:)), hold on;
     plot(x,ana(floor(0.5*t_res),:)), hold on;
     plot(x,-td*x), hold off;
     axis([-0.2 6 -0.03 0.04]);
     title(sprintf('Run-down ($t=%0.2f$)',Th), 'Interpreter','latex');
     xlabel('$x$', 'Interpreter','latex');
     ylabel('$\eta$', 'Interpreter','latex');
     legend('Numerical', 'Analytical','Bathymetry');

     subplot(4,1,4);
     plot(t,stat_norm);
     title('L2 Norm', 'Interpreter','latex');
     xlabel('$t$', 'Interpreter','latex');
     ylabel('L2 Norm', 'Interpreter','latex');
     cnt = cnt+1;
   else
     ;
   end

   % Wave Animation
   if contains(Wave_Animation,"on")
     figure(cnt)
     for i = 1:t_res
       plot(x,num(i,:) ), hold on;
       plot(x,ana(i,:) ), hold on;
       plot(x,-td*x), hold off;
       axis([x0 Xf -0.05 0.05])
       pause(0.05);
     end
   else
     ;
   end
  plots = 1;
end
