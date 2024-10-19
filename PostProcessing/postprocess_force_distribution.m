% PLOT_FORCE_DISTRIBUTION
% This script plots the distribution of forces in the x and y directions from F_store at 5 time points.
%
% Variables:
%   F_store - 3D array storing force data [time_points, Npts_t, 2] where:
%             - 3rd dimension = 1 is force in the x-direction (F_x)
%             - 3rd dimension = 2 is force in the y-direction (F_y)

% Check if the required variable F_store is available in the workspace
if ~exist('F_store', 'var')
    error('F_store variable is missing.');
end

% Determine the dimensions of F_store
[time_points, Npts_t, ~] = size(F_store);

% Select 5 time points evenly spaced throughout the simulation
selected_time_steps = round(linspace(1, time_points, 5));

% Number of bins for the histograms
num_bins = 15;  % Adjust this value as needed

% Loop through the selected time steps to plot the force distribution
for k = 1:length(selected_time_steps)
    t = selected_time_steps(k);
    
    % Extract the forces at time step t
    F_x = F_store(t, :, 1); % Forces in x-direction
    F_y = F_store(t, :, 2); % Forces in y-direction
    
    % Create a new figure for the current time step
    figure('Name', ['Force Distribution at Time Step ', num2str(t)], 'NumberTitle', 'off');
    
    % Subplot 1: Histogram of forces in the x-direction
    subplot(2, 1, 1);
    histogram(F_x, num_bins, 'Normalization', 'pdf'); % PDF normalization with increased bins
    xlabel('Force in X direction (N)');
    ylabel('Probability Density');
    title(['Force Distribution (X) at Time Step ', num2str(t)]);
    grid on;
    
    % Subplot 2: Histogram of forces in the y-direction
    subplot(2, 1, 2);
    histogram(F_y, num_bins, 'Normalization', 'pdf'); % PDF normalization with increased bins
    xlabel('Force in Y direction (N)');
    ylabel('Probability Density');
    title(['Force Distribution (Y) at Time Step ', num2str(t)]);
    grid on;
end

disp('Force distribution plots for selected time steps generated successfully.');
