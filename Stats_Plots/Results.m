function [runup_ana, runup_num] = Results(stat_norm_max,stat_norm_max_t);

  global runup_num_i runup_ana_i stat_norm ana num
  global x x0 Xf x_res t Tf t0 t_res display_stat
  global runup_ana_t runup_num_t runup_ana_max_i runup_num_max_i


  for i = 1:t_res
    runup_ana(:,i) = ana(i,runup_ana_i(:,i));    % eta at wet/dry boundary
    runup_num(:,i) = num(i,runup_num_i(:,i));
  end

  runup_ana_max = max(runup_ana);                 % maximum value of eta
  runup_num_max = max(runup_num);
  runup_ana_max_i = find(runup_ana == max(runup_ana));
  runup_num_max_i = find(runup_num == max(runup_num));
  runup_ana_t = runup_ana_max_i*((Tf-t0)/(t_res-1));      % time of maximum runup
  runup_num_t = runup_num_max_i*((Tf-t0)/(t_res-1));
  runup_ana_x = x0+max(runup_ana_i)*((Xf-x0)/(x_res-1));  % position of maximum runup
  runup_num_x = x0+max(runup_num_i)*((Xf-x0)/(x_res-1));

  % displaying results
  fprintf(repmat('\b',1,display_stat));
  fprintf('\n');
  fprintf('Max L2 norm: %0.6f at t=%0.3f \n',stat_norm_max,...
    stat_norm_max_t)
  fprintf('Max Run-up of Analytical Solution: eta = %0.4f at (x,t) = (%0.4f,%0.4f).\n',...
    runup_ana_max,runup_ana_x,runup_ana_t);
  fprintf('Max Run-up of Numerical Solution: eta = %0.4f at (x,t) = (%0.4f,%0.4f).\n',...
    runup_num_max,runup_num_x,runup_num_t);

end
