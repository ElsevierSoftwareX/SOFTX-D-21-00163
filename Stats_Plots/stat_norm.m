function [ana, ana_u, num, num_u]= stat_norm(eta_analytic, eta_fvm, u_analytic, u_fvm)

  % global variables
  global x x0 Xf x_res t t_res g l td display_stat
  global num ana num_u ana_u runup_num_i runup_ana_i stat_norm

  fprintf('\n');
  display_stat = fprintf('Post-processing solution arrays.');

  % initializing processsed solution arrays
  [ana, ana_u, num, num_u] = deal(zeros(t_res,x_res));
  [runup_ana_i, runup_num_i] = deal(zeros(1,t_res));

  for i = 1:t_res
    % assigning values to empty arrays
    ana(i,:) = eta_analytic(x, repmat(t(i), 1, x_res)).*(l*td);
    ana_u(i,:) = u_analytic(x, repmat(t(i), 1, x_res)).*sqrt((g*td)/l);
    num(i,:) = eta_fvm(:,i)'.*(l*td);
    num_u(i,:) = u_fvm(:,i)'.*sqrt((g*td)/l);
    maximum = 0;
    for j = 1:x_res
      if ana(i, j) + td*x(j) < 10^(-6)   % processing eta (analytical)
        ana(i, j) = NaN;
        maximum_ana = j;  % index of wet/dry boundary
        runup_ana_i(1,i) = maximum_ana+1;
      end
      if num(i, j) + td*x(j) < 0   % processing eta (numerical)
        num(i, j) = NaN;
        maximum_num = j;  % index of wet/dry boundary
        runup_num_i(1,i) = maximum_num+1;
      end
      if ana_u(i, j) + td*x(j) < 0   % processing u (analytical)
        ana_u(i, j) = NaN;
      end
      if num_u(i, j) + td*x(j) < 0   % processing u (numerical)
        num_u(i, j) = NaN;
      end
    end
    if i == 1
        normalize = norm(ana(1, maximum_ana+1:end));
    end
    if maximum_ana > maximum_num  % Computing l2 norm for all t
      stat_norm(i) = (norm(ana(i, maximum_ana+1:end) - num(i, maximum_ana+1:end))...
        *(Xf-x0)/(x_res-1))/(normalize...
        *(Xf-x0)/(x_res-1));
    elseif maximum_num > maximum_ana
      stat_norm(i) = (norm(ana(i, maximum_num+1:end) - num(i, maximum_num+1:end))...
        *(Xf-x0)/(x_res-1))/(normalize...
        *(Xf-x0)/(x_res-1));
    elseif maximum_num == maximum_ana
      stat_norm(i) = (norm(ana(i, maximum_num+1:end) - num(i, maximum_num+1:end))...
        *(Xf-x0)/(x_res-1))/(normalize...
        *(Xf-x0)/(x_res-1));
    end
  end
end
