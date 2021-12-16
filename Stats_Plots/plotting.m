function plots = plotting(WaveAnimation, Eta_3D, U_3D, Runup, L2Norm, InitialConditions,...
  DataProjection)

  % function that produces all plots
  % variables needed for plotting
  global x x0 Xf t Tf t_res td n
  global proj0 proj1 projn s u0 eta0
  global num num_u ana ana_u stat_norm
  global runup_ana_t runup_num_t runup_ana_max_i runup_num_max_i

  stat_norm_max_i = find(stat_norm == max(stat_norm));   % index with max L2 norm
  stat_norm_max_t = Tf/t_res*stat_norm_max_i;            % time with max l2 norm
  stat_norm_max = stat_norm(stat_norm_max_i);            % max L2 norm

  [runup_ana, runup_num] = Results(stat_norm_max,stat_norm_max_t);

  cnt = 1;   % figure counter

  % Wave Animation
  if contains(WaveAnimation,"on")
    figure(cnt)
    for i = 1:t_res
      plot(x,num(i,:) ), hold on;
      plot(x,ana(i,:) ), hold on;
      plot(x,-td*x,'k','LineWidth',2), hold off;
      axis([x0 Xf -0.05 0.05])
      title('Runup','Interpreter','latex','Fontsize', 14)
      xlabel('$x$','Interpreter','latex','Fontsize', 14)
      ylabel('$\eta$','Interpreter','latex','Fontsize', 14)
      legend({'Numerical','Analytical'},'Location','southwest','Interpreter','latex');
      drawnow limitrate
    end
    drawnow
    cnt = cnt+1;
  else
    ;
  end

  % eta(x,t)
  if contains(Eta_3D,"on")
    figure(cnt)
    surf(x,t,num)
    colormap winter
    shading interp
    title('Wave Height (Numerical)','Interpreter','latex','Fontsize', 14)
    xlabel('$x$','Interpreter','latex','Fontsize', 14)
    ylabel('$t$','Interpreter','latex','Fontsize', 14)
    zlabel('$\eta$','Interpreter','latex','Fontsize', 14)
    cnt = cnt+1;

    figure(cnt)
    surf(x,t,ana)
    colormap winter
    shading interp
    title('Wave Height (Analytical)','Interpreter','latex','Fontsize', 14)
    xlabel('$x$','Interpreter','latex','Fontsize', 14)
    ylabel('$t$','Interpreter','latex','Fontsize', 14)
    zlabel('$\eta$','Interpreter','latex','Fontsize', 14)
    cnt = cnt+1;
  else
    ;
  end

  % u(x,t)
  if contains(U_3D,"on")
    figure(cnt)
    surf(x,t,num_u)
    colormap parula
    shading interp
    title('Wave Velocity (Numerical)','Interpreter','latex','Fontsize', 14)
    xlabel('$x$','Interpreter','latex','Fontsize', 14)
    ylabel('$t$','Interpreter','latex','Fontsize', 14)
    zlabel('$u$','Interpreter','latex','Fontsize', 14)
    cnt = cnt+1;

    figure(cnt)
    surf(x,t,ana_u)
    colormap parula
    shading interp
    title('Wave Velocity (Analytical)','Interpreter','latex','Fontsize', 14)
    xlabel('$x$','Interpreter','latex','Fontsize', 14)
    ylabel('$t$','Interpreter','latex','Fontsize', 14)
    zlabel('$u$','Interpreter','latex','Fontsize', 14)
    cnt = cnt+1;
  else
    ;
  end

  % Data Projection
  if contains(DataProjection,"on")
    if n == 1
      figure(cnt);
      plot(s,proj1(1,:),'k'), hold on;
      plot(s, u0(x),'--k'), hold off;
      title('Projection of $\varphi_0$', 'Interpreter','latex','Fontsize', 14);
      xlabel('$\sigma$', 'Interpreter','latex','Fontsize', 14);
      ylabel('$\varphi$', 'Interpreter','latex','Fontsize', 14);
      legend({'$1^{st}$ order projection',...
        '$\varphi_0$'},'Location','southwest','Interpreter','latex');
        cnt = cnt+1;

        figure(cnt)
        plot(s,proj1(2,:)), hold on;
        plot(s, eta0(s)+u0(s).^2/2), hold off;
        title('Projection of $\psi_0$', 'Interpreter','latex','Fontsize', 14);
        xlabel('$\sigma$', 'Interpreter','latex','Fontsize', 14);
        ylabel('$\psi_n$', 'Interpreter','latex','Fontsize', 14);
        legend({'$1^{st}$ order projection',...
        '$\psi_0$'},'Location','southwest','Interpreter','latex');
        cnt = cnt+1;
    else
      figure(cnt);
      plot(s,projn(1,:),'k'), hold on;
      plot(s,proj1(1,:),'--k'), hold on;
      plot(s, u0(s),':k'), hold off;
      title('Projection of $\varphi_0$', 'Interpreter','latex','Fontsize', 14);
      xlabel('$\sigma$', 'Interpreter','latex','Fontsize', 14);
      ylabel('$\varphi$', 'Interpreter','latex','Fontsize', 14);
      legend({sprintf('%dth order projection',n)','$1^{st}$ order projection',...
        '$\varphi_0$'},'Location','southwest','Interpreter','latex');

      cnt = cnt+1;
      figure(cnt)

      plot(s,projn(2,:),'k'), hold on;
      plot(s,proj1(2,:),'--k'), hold on;
      plot(s, eta0(2)+(u0(2).^2)/2,':k'), hold off;
      title('Projection of $\psi_0$', 'Interpreter','latex','Fontsize', 14);
      xlabel('$\sigma$', 'Interpreter','latex','Fontsize', 14);
      ylabel('$\psi_n$', 'Interpreter','latex','Fontsize', 14);
      legend({sprintf('%dth order projection', n),'$1^{st}$ order projection',...
        '$\psi_0$'},'Location','southwest','Interpreter','latex');

      cnt = cnt+1;
    end
   else
     ;
   end

   if contains(Runup,"on")
     figure(cnt)
     plot(t,runup_ana,'k'), hold on;
     plot(t,runup_num,':k'), hold off;
     title('Maximum Runup Values', 'Interpreter','latex','Fontsize', 14);
     xlabel('$t$', 'Interpreter','latex','Fontsize', 14);
     ylabel('$\eta$', 'Interpreter','latex','Fontsize', 14);
     legend({'Analytical','Numerical'},'Location','southwest','Interpreter','latex');
     cnt = cnt+1;
   else
     ;
   end

   % L2 norm plot
   if contains(L2Norm,"on")

     figure(cnt)
     subplot(4,1,1);
     plot(x,num(runup_num_max_i,:),'--k'), hold on;
     plot(x,ana(runup_ana_max_i,:),'k'), hold on;
     plot(x,-td*x,'k','LineWidth',2), hold off;
     axis([-0.2 6 -0.03 0.04]);
     title(sprintf('Maximum Run-up (Analytical: $t=%0.2f$ Numerical: $t=%0.2f$)',...
     runup_ana_t,runup_num_t), 'Interpreter','latex','Fontsize', 14);
     xlabel('$x$', 'Interpreter','latex');
     ylabel('$\eta$', 'Interpreter','latex');
     legend('Numerical', 'Analytical','Bathymetry');

     subplot(4,1,2);
     plot(x,num(stat_norm_max_i,:),'--k'), hold on;
     plot(x,ana(stat_norm_max_i,:),'k'), hold on;
     plot(x,-td*x,'k','LineWidth',2), hold off;
     axis([-0.2 6 -0.03 0.04]);
     title(append('$\eta$ ', sprintf('at Time of Maximum L2 Norm ($t=%0.2f$)',stat_norm_max_t)),...
      'Interpreter','latex','Fontsize', 14);
     xlabel('$x$', 'Interpreter','latex');
     ylabel('$\eta$', 'Interpreter','latex');
     legend('Numerical', 'Analytical','Bathymetry');

     subplot(4,1,3);
     plot(x,num_u(stat_norm_max_i,:),'--k'), hold on;
     plot(x,ana_u(stat_norm_max_i,:),'k'), hold on;
     plot(x,-td*x,'k','LineWidth',2), hold off;
     axis([-0.2 6 -0.03 0.04]);
     title(sprintf('$u$ at Time of Maximum L2 Norm ($t=%0.2f$)',stat_norm_max_t),...
      'Interpreter','latex','Fontsize', 14);
     xlabel('$x$', 'Interpreter','latex');
     ylabel('$\eta$', 'Interpreter','latex');
     legend('Numerical', 'Analytical','Bathymetry');

     subplot(4,1,4);
     plot(t,stat_norm,'k');
     title('L2 Norm (normalized)', 'Interpreter','latex','Fontsize', 14);
     xlabel('$t$', 'Interpreter','latex','Fontsize', 14);
     ylabel('L2 Norm', 'Interpreter','latex','Fontsize', 14);
     cnt = cnt+1;
   else
     ;
   end

   % initial conditions
   if contains(InitialConditions,"on")
     figure(cnt)
     plot(x,num(1,:),'k')
     title('Initial Wave Profile ($\eta$)', 'Interpreter','latex','Fontsize', 14);
     xlabel('$x$', 'Interpreter','latex','Fontsize', 14);
     ylabel('$t$', 'Interpreter','latex','Fontsize', 14);
     cnt = cnt+1;

     figure(cnt)
     plot(x,num_u(1,:),'k')
     title('Initial Wave Velocity ($u$)', 'Interpreter','latex','Fontsize', 14);
     xlabel('$x$', 'Interpreter','latex','Fontsize', 14);
     ylabel('$t$', 'Interpreter','latex','Fontsize', 14);
     cnt = cnt+1;
   end

%}
  plots = cnt;
end
