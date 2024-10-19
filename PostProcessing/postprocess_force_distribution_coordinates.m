% PLOT_FORCE_VECTORS_WITH_VALUES
% This script plots the force vectors using the coordinates in mat_r and the forces in F_store,
% and annotates the plot with the magnitudes or components of the forces.
%
% Variables:
%   F_store - 3D array storing force data [time_points, Npts_t, 2] where:
%             - 3rd dimension = 1 is force in the x-direction (F_x)
%             - 3rd dimension = 2 is force in the y-direction (F_y)
%   mat_r   - 2D array storing the coordinates of the points [Npts_t, 2]
%             - 1st column = x-coordinates
%             - 2nd column = y-coordinates

% Check if the required variables F_store and mat_r are available in the workspace
if ~exist('F_store', 'var') || ~exist('mat_r', 'var')
    error('F_store and/or mat_r variables are missing.');
end

% Determine the dimensions of F_store
[time_points, Npts_t, ~] = size(F_store);

% Select a specific time step for the plot (e.g., the last time step)
t = time_points;  % You can change this to another time step if desired

% Extract the forces at the selected time step
F_x = F_store(t, :, 1); % Forces in x-direction
F_y = F_store(t, :, 2); % Forces in y-direction

% Extract the coordinates from mat_r
x_coords = mat_r(:, 1)'; % x-coordinates
y_coords = mat_r(:, 2)'; % y-coordinates

% Calculate the magnitude of the force vectors
F_magnitude = sqrt(F_x.^2 + F_y.^2);

% Calculate the ratio |F_x| / |F_y| for each vector (handling F_y = 0)
F_ratio = abs(F_x) ./ max(abs(F_y), eps); % eps prevents division by zero

% Find the index of the vector with the maximum ratio and largest magnitude
[~, most_horizontal_idx] = max(F_ratio .* F_magnitude); % Combined criterion

% Get the most horizontal and largest force vector
most_horizontal_force_x = F_x(most_horizontal_idx);
most_horizontal_force_y = F_y(most_horizontal_idx);
most_horizontal_magnitude = F_magnitude(most_horizontal_idx);

% Get the coordinates of this force vector
x_coord = mat_r(most_horizontal_idx, 1);
y_coord = mat_r(most_horizontal_idx, 2);

% Display the result
fprintf('Most Horizontal and Largest Force Vector:\n');
fprintf('X Coordinate: %.2f\n', x_coord);
fprintf('Y Coordinate: %.2f\n', y_coord);
fprintf('Force in X Direction (F_x): %.2f N\n', most_horizontal_force_x);
fprintf('Force in Y Direction (F_y): %.2f N\n', most_horizontal_force_y);
fprintf('Force Magnitude: %.2f N\n', most_horizontal_magnitude);

% Create a new figure for the quiver plot with white background
figure('Name', ['Force Vector Plot with Highlighted Vector at Time Step ', num2str(t)], ...
       'NumberTitle', 'off', 'Color', 'white'); % Set background color to white

% Generate the quiver plot for all forces
quiver(x_coords, y_coords, F_x, F_y, 'AutoScale', 'on', 'MaxHeadSize', 0.5);

% Highlight the most horizontal and largest vector in red
hold on;
quiver(x_coord, y_coord, most_horizontal_force_x, most_horizontal_force_y, ...
       'AutoScale', 'on', 'MaxHeadSize', 0.5, 'Color', 'r', 'LineWidth', 2);

% Add labels and title
xlabel('X Coordinate', 'FontSize', 14);
ylabel('Y Coordinate', 'FontSize', 14);
title(['Force Vector Plot with Highlighted Vector at Time Step ', num2str(t)], 'FontSize', 14);
grid on;

% Optionally, adjust the axis to make the vectors more visible
axis equal;  % Ensure the axes have the same scale
xlim([min(x_coords), max(x_coords)]);
ylim([min(y_coords), max(y_coords)]);

% Annotate the vectors with the magnitude of the force (font size 14)
% for i = 1:Npts_t
%     % Annotate with magnitude (force value)
%     text(x_coords(i), y_coords(i), ['|F|=', num2str(F_magnitude(i), '%.2f')], ...
%         'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
% end

% Annotate the highlighted vector with its magnitude
text(x_coord, y_coord, ['|F|=', num2str(most_horizontal_magnitude, '%.2f')], ...
    'FontSize', 14, 'Color', 'r', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');

hold off;

disp('Force vector plot with highlighted vector generated successfully.');