function stats = stat_norm()

  % global variables
  global x x0 Xf x_res t t_res Tf t0 td num ana eta_analytic eta_fvm
  global stat_norm_max_i stat_norm_max_t stat_norm
  global max_num_t max_ana_t max_num_t_i max_ana_t_i t_offset max_ana_x

  % Initialing solution arrays
  num = zeros(t_res, x_res);
  ana = zeros(t_res, x_res);

  t_offset = 5;

  fprintf('\n');
  display_stat = fprintf('Post-processing solution arrays.');
  for i = t_offset:t_res

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
    % Computing L2 norm for all t
    stat_norm(i) = norm(ana(i, maximum+1:end) - num(i, maximum+1:end));
  end

  fprintf(repmat('\b',1,display_stat));
  display_stat = fprintf('Mean and Maximum L2 norm.');

  stat_norm_max_i = find(stat_norm == max(stat_norm));   % index with max L2 norm
  stat_norm_max_t = Tf/t_res*stat_norm_max_i;            % time with max l2 norm
  stat_norm_max = stat_norm(stat_norm_max_i);            % max L2 norm
  stat_norm_avg = mean(stat_norm(:));                    % mean L2 norm

  fprintf(repmat('\b',1,display_stat));
  display_stat = fprintf('Processing run-up data.');

  [max_num_t_i, max_num_x_i] = find(num == max(num,[],'all'));   % Finding maximum run-up in solution arrays
  [max_ana_t_i, max_ana_x_i] = find(ana == max(ana,[],'all'));   % ...
  max_num_t = (max_num_t_i+t_offset-1)*((Tf-t0)/(t_res-1));      % Produces time of maximum runup
  max_ana_t = (max_ana_t_i+t_offset-1)*((Tf-t0)/(t_res-1));      % ...
  max_num_x = x0+(max_num_x_i-1)*((Xf-x0)/(x_res-1));            % Produces position of maximum runup
  max_ana_x = x0+(max_ana_x_i-1)*((Xf-x0)/(x_res-1));            % ...
  x_diff = abs((max_ana_x-max_num_x)/((max_ana_x+max_num_x)/2))*100;   % Percent difference in two solutions
  t_diff = abs((max_ana_t-max_num_t)/((max_ana_t+max_num_t)/2))*100;   % ...
  x_error = (Xf-x0)/(x_res-1);
  t_error = ((Tf-t0)/(t_res-1))*100;

  if max_num_x > max_ana_x                    % Error calculations
    x_error = abs((x_error/max_ana_x))*100;
  else
    x_error = abs((x_error/max_num_x))*100;
  end

  % Displaying results
  fprintf(repmat('\b',1,display_stat));
  fprintf('\n');
  fprintf('Max L2 norm: %0.6f at t=%0.3f \n', stat_norm_max,...
    stat_norm_max_t);
  fprintf('Mean L2 norm: %0.6f \n', stat_norm_avg);
  fprintf('\n');
  fprintf('Max Run-up of Analytical Solution: x=%0.6f at t=%0.6f. \n', max_ana_x,...
    max_ana_t);
  fprintf('Max Run-up of Numerical Solution: x=%0.4f at t=%0.3f. \n', max_num_x,...
    max_num_t);
  fprintf('\n');
  fprintf('Difference in Position of Maximum Run-up: %0.3f%%. \n', x_diff);
  fprintf('Difference in Time of Maximum Run-up: %0.3f%%. \n', t_diff);
  fprintf('\n');
  fprintf('Error: x<+%0.3f%% and t=±%0.03f%%. \n', x_error, t_error);
  stats=1;

end
