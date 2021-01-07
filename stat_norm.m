function stats = stat_norm()

  % global variables
  global x x_res t t_res Tf td num ana eta_analytic eta_fvm
  global stat_norm_max_i stat_norm_max_t stat_norm

  num = zeros(t_res, x_res);
  ana = zeros(t_res, x_res);

  for i = 1:t_res

    ana(i,:) = eta_analytic(x, repmat(t(i), 1, x_res));
    num(i,:) = eta_fvm(:, i)';

    % keeping wave above beach in plotting
    maximum = 0;
    for j = 1:x_res
      if num(i, j) + td*x(j) < 0
        num(i, j) = NaN;
        maximum = j;
      end
      if ana(i,j) + td*x(j) < 0
        ana(i,j) = NaN;
        maximum = j;
      end
    end
    stat_norm(i) = norm(ana(i, maximum+1:end) - num(i, maximum+1:end));
  end

  stat_norm_max_i = find(stat_norm == max(stat_norm));   % index with max L2 norm
  stat_norm_max_t = Tf/t_res*stat_norm_max_i;            % time with max l2 norm
  stat_norm_max = stat_norm(stat_norm_max_i);            % max L2 norm
  stat_norm_avg = mean(stat_norm(:));                    % mean L2 norm

  fprintf('\n');
  fprintf('Max L2 norm: %0.6f at t = %0.2f \n', stat_norm_max,...
   stat_norm_max_t);
  fprintf('Mean L2 norm: %0.6f \n', stat_norm_avg);
  stats=1;
end
