% PLOT_SIMULATION_RESULTS
% This script takes the stored simulation results and plots them over time.
%
% Variables:
%   F_t_res              - Cell array of force validation results over time
%   F_adhesion_res       - Cell array of adhesion forces over time
%   force_difference_res - Cell array of force differences over time
%   force_ratio_res      - Cell array of force ratios over time
%   time_store           - (Optional) Array of time points

% Check if the required variables are available in the workspace
if ~exist('F_t_res', 'var') || ~exist('F_adhesion_res', 'var') || ...
    ~exist('force_difference_res', 'var') || ~exist('force_ratio_res', 'var')
error('One or more required variables (F_t_res, F_adhesion_res, force_difference_res, force_ratio_res) are missing.');
end

% Determine the dimensions of the input data
[nRows, nCols] = size(F_t_res);

% Determine the length of time points from one of the result arrays
nTimePoints = length(F_t_res{1,1});

% If time_store does not exist, generate a default time array (1, 2, 3, ..., nTimePoints)
if ~exist('time_store', 'var')
disp('time_store not found. Generating default time array...');
time_store = 1:nTimePoints; % Default time array
end

% Check if F_t_res contains non-empty data
if any(~cellfun(@isempty, F_t_res))
figure('Name', 'Force Validation (F_t)', 'NumberTitle', 'off');
hold on;
end

if any(~cellfun(@isempty, F_adhesion_res))
figure('Name', 'Adhesion Force (F_adhesion)', 'NumberTitle', 'off');
hold on;
end

if any(~cellfun(@isempty, force_difference_res))
figure('Name', 'Force Difference', 'NumberTitle', 'off');
hold on;
end

if any(~cellfun(@isempty, force_ratio_res))
figure('Name', 'Force Ratio', 'NumberTitle', 'off');
hold on;
end

% Combined plot for F_t, F_adhesion, and force_difference
combined_plot_needed = any(~cellfun(@isempty, F_t_res)) && any(~cellfun(@isempty, F_adhesion_res)) && any(~cellfun(@isempty, force_difference_res));
if combined_plot_needed
figure('Name', 'Combined Plot (F_t, F_adhesion, Force Difference)', 'NumberTitle', 'off');
hold on;
end

% Loop through each simulation result and plot the data
for i = 1:nRows
for j = 1:nCols
    % Plot Force Validation (F_t)
    if ~isempty(F_t_res{i,j})
        figure(1); % Switch to Force Validation plot
        plot(time_store, F_t_res{i,j}, 'DisplayName', ['Sim(', num2str(i), ',', num2str(j), ')']);
        xlabel('Time (s)');
        ylabel('Force Validation (N)');
        title('Force Validation (F_t) over Time');
        legend show;
    end
    
    % Plot Adhesion Force (F_adhesion)
    if ~isempty(F_adhesion_res{i,j})
        figure(2); % Switch to Adhesion Force plot
        plot(time_store, F_adhesion_res{i,j}, 'DisplayName', ['Sim(', num2str(i), ',', num2str(j), ')']);
        xlabel('Time (s)');
        ylabel('Adhesion Force (N)');
        title('Adhesion Force (F_adhesion) over Time');
        legend show;
    end

    % Plot Force Difference
    if ~isempty(force_difference_res{i,j})
        figure(3); % Switch to Force Difference plot
        plot(time_store, force_difference_res{i,j}, 'DisplayName', ['Sim(', num2str(i), ',', num2str(j), ')']);
        xlabel('Time (s)');
        ylabel('Force Difference (N)');
        title('Force Difference over Time');
        legend show;
    end

    % Plot Force Ratio
    if ~isempty(force_ratio_res{i,j})
        figure(4); % Switch to Force Ratio plot
        plot(time_store, force_ratio_res{i,j}, 'DisplayName', ['Sim(', num2str(i), ',', num2str(j), ')']);
        xlabel('Time (s)');
        ylabel('Force Ratio');
        title('Force Ratio over Time');
        legend show;
    end

    % Plot F_t, F_adhesion, and force_difference together
    if combined_plot_needed && ~isempty(F_t_res{i,j}) && ~isempty(F_adhesion_res{i,j}) && ~isempty(force_difference_res{i,j})
        figure(5); % Combined plot
        plot(time_store, F_t_res{i,j}, 'r', 'DisplayName', ['F_t (Sim ', num2str(i), ',', num2str(j), ')']);
        hold on;
        plot(time_store, F_adhesion_res{i,j}, 'g', 'DisplayName', ['F_adhesion (Sim ', num2str(i), ',', num2str(j), ')']);
        plot(time_store, force_difference_res{i,j}, 'b', 'DisplayName', ['Force Difference (Sim ', num2str(i), ',', num2str(j), ')']);
        xlabel('Time (s)');
        ylabel('Force (N)');
        title('Combined Plot (F_t, F_adhesion, Force Difference)');
        legend show;
    end
end
break;
end

% Hold off after plotting (for all figures that were created)
if any(~cellfun(@isempty, F_t_res))
hold off;
end

if any(~cellfun(@isempty, F_adhesion_res))
hold off;
end

if any(~cellfun(@isempty, force_difference_res))
hold off;
end

if any(~cellfun(@isempty, force_ratio_res))
hold off;
end

if combined_plot_needed
hold off;
end

disp('Plots generated successfully.');
