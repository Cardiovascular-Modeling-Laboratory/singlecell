% run_simulations_with_OOP_calculation.m
% Run simulations with multiple max_stretch_factor and fold_time,
% and calculate the Orientation Order Parameter.
% If an error occurs, save the error log and continue to the next simulation.

% Setting max_stretch_factors and fold_time
% max_stretch_factors = [3.5, 4.0, 5.0]; % List of max_stretch_factors to use
% fold_times = [(72-1)*3600, (72-3)*3600, (72-6)*3600, (72-12)*3600, (72-24)*3600]; % List of fold_times to use
% fold_times = [(72-1)*3600]; % List of fold_times to use
max_stretch_factors = [1.0, 2.0, 3.0, 4.0, 5.0]; % List of max_stretch_factors to use
nuc_rel_vec=[0, 1];
n_trial = 5; % Number of trials for each condition

% Directory to save error logs
error_log_dir = fullfile(pwd, 'error_logs');
if ~exist(error_log_dir, 'dir')
    mkdir(error_log_dir);
end

% Loop through each max_stretch_factor and run simulations
for i = 1:length(max_stretch_factors)
    max_stretch_factor = max_stretch_factors(i);
    for j = 1:length(nuc_rel_vec)
        nuc_rel = nuc_rel_vec(j);

        for k = 1:n_trial
            try
                % Call single_cell_units_linked_v3 and save the results
                disp(['Running simulation with max_stretch_factor = ', num2str(max_stretch_factor), ...
                    ' and nuc_rel = ', num2str(nuc_rel)]);
                single_cell_units_linked_v3_NucPlacement(max_stretch_factor, [nuc_rel], k);

                % Calculate the Orientation Order Parameter
                % calculate_orientation_order_parameter;
                
            catch ME
                % Save error information
                errorTime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
                errorMessage = ME.message;
                errorDetails = struct('Time', errorTime, 'Message', errorMessage, 'Stack', ME.stack);
                error_filename = fullfile(error_log_dir, ['error_', errorTime, '_factor_', ...
                                num2str(max_stretch_factor), '_fold_', num2str(nuc_rel), '.mat']);
                save(error_filename, 'errorDetails');
                disp(['Error occurred and saved to ', error_filename]);
                disp(['Error message: ', errorMessage]);
                disp(['Error stack: ']);
                ME.stack
            end
        end
    end
    close all;
end
