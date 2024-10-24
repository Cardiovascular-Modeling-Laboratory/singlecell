function postprocess_force_distribution_multisim()
    % COMPARE_FORCE_DISTRIBUTION
    % This function compares the force distributions in the x- and y-directions
    % across multiple simulation results stored in .mat files.
    % Files are selected manually through a file dialog.
    %
    % Usage:
    %   compare_force_distribution()
    %   The function will open a file dialog to select .mat files for comparison.
    
    % Open file selection dialog and allow selection of multiple files
    [file_names, file_path] = uigetfile('*.mat', 'Select One or More MAT Files', 'MultiSelect', 'on');
    
    % If no file is selected, exit the function
    if isequal(file_names, 0)
        disp('No files were selected. Exiting.');
        return;
    end
    
    % Convert file_names to a cell array if only one file is selected
    if ischar(file_names)
        file_names = {file_names};
    end
    
    % Number of files (simulations) to compare
    num_files = length(file_names);

    % Prepare colors for distinguishing the plots from different simulations
    colors = lines(num_files);  % 'lines' generates distinct colors

    % Create a figure for force distribution comparison
    figure('Name', 'Force Distribution Comparison', 'NumberTitle', 'off');
    
    % Initialize subplots for x- and y-direction forces
    subplot(2, 1, 1);  % Subplot for x-direction forces
    hold on;
    xlabel('Force in X direction (N)');
    ylabel('Probability Density');
    title('Force Distribution in X Direction');
    grid on;
    
    subplot(2, 1, 2);  % Subplot for y-direction forces
    hold on;
    xlabel('Force in Y direction (N)');
    ylabel('Probability Density');
    title('Force Distribution in Y Direction');
    grid on;

    % Iterate through each selected file, load the data, and plot the force distributions
    for k = 1:num_files
        filename = fullfile(file_path, file_names{k});

        % Load the .mat file (assume it contains 'F_store' and 'mat_r')
        data = load(filename, 'F_store', 'mat_r');
        
        % Check if F_store exists in the file
        if ~isfield(data, 'F_store')
            error('The file %s does not contain F_store.', filename);
        end
        
        % Get the dimensions of F_store
        [time_points, Npts_t, ~] = size(data.F_store);

        % Select the last time step for comparison (you can change this)
        t = time_points;  % Last time step
        
        % Extract forces at the selected time step
        F_x = data.F_store(t, :, 1);  % Forces in x-direction
        F_y = data.F_store(t, :, 2);  % Forces in y-direction
        
        % Plot the force distributions in x-direction (subplot 1)
        subplot(2, 1, 1);
        histogram(F_x, 50, 'Normalization', 'pdf', 'DisplayName', file_names{k}, ...
            'EdgeColor', colors(k, :), 'LineWidth', 1.5);

        % Plot the force distributions in y-direction (subplot 2)
        subplot(2, 1, 2);
        histogram(F_y, 50, 'Normalization', 'pdf', 'DisplayName', file_names{k}, ...
            'EdgeColor', colors(k, :), 'LineWidth', 1.5);
    end

    % Add legends to each subplot
    subplot(2, 1, 1);
    legend show;

    subplot(2, 1, 2);
    legend show;

    hold off;
end
