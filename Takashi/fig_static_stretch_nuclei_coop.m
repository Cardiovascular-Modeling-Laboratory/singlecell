% This script calculates COOP between datasets with and without nuclei and plots the results as a bar graph.

%% Set data directory and file pattern
dataDir = pwd;
filePattern = fullfile(dataDir, 'Takashi', 'results', 'fig9_v2_n=*_final*_nuc_rel_vec_*.mat');
resultFiles = dir(filePattern);

% Initialize file information
file_info = struct('filename', {}, 'n_sim', {}, 'final_value', {}, 'nucleus_obstruction', {});

% Collect file information
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    % Extract necessary information from the filename
    tokens = regexp(resultFiles(k).name, 'fig9_v2_n=(\d+)_final(\d+)_nuc_rel_vec_(\d+)', 'tokens');
    if ~isempty(tokens)
        n_sim = str2double(tokens{1}{1});
        final_value = str2double(tokens{1}{2});
        nucleus_obstruction = str2double(tokens{1}{3});
        file_info(end+1) = struct('filename', filename, 'n_sim', n_sim, 'final_value', final_value, 'nucleus_obstruction', nucleus_obstruction);
    else
        disp(['Could not parse filename: ', resultFiles(k).name]);
        continue;
    end
end

% Get unique stretch conditions
final_values = [file_info.final_value];
unique_final_values = unique(final_values);

% Initialize COOP results
COOP_results = [];

% Process each stretch condition
for i = 1:length(unique_final_values)
    fv = unique_final_values(i);
    % Get files with and without nuclei
    idx_no_nucleus = find([file_info.final_value] == fv & [file_info.nucleus_obstruction] == 0);
    idx_with_nucleus = find([file_info.final_value] == fv & [file_info.nucleus_obstruction] == 1);

    % Determine the number of pairs
    num_pairs = min(length(idx_no_nucleus), length(idx_with_nucleus));
    for j = 1:num_pairs
        % Load data without nuclei
        file_no = file_info(idx_no_nucleus(j)).filename;
        [OOP_no, grid_no] = calculate_OOP_and_fiber_orientation_grid(file_no);

        % Load data with nuclei
        file_with = file_info(idx_with_nucleus(j)).filename;
        [OOP_with, grid_with] = calculate_OOP_and_fiber_orientation_grid(file_with);

        % Calculate COOP
        COOP = calculate_COOP(OOP_no, OOP_with, grid_no, grid_with);
        % Record the presence or absence of nuclei in the results
        COOP_results = [COOP_results; struct('final_value', fv, 'nucleus_obstruction', 0, 'COOP', COOP)];
    end
end

% Convert results to data table
COOP_table = struct2table(COOP_results);

% Sort data
[sorted_final_values, sort_idx] = sort(COOP_table.final_value);
sorted_COOP = COOP_table.COOP(sort_idx);
sorted_nucleus_obstruction = COOP_table.nucleus_obstruction(sort_idx);

% Convert group variable to categorical
FinalValue = categorical(sorted_final_values);

% Statistical test by ANOVA and Tukey HSD (excluding NucleusObstruction)
group = FinalValue;
[p, tbl, stats] = anovan(sorted_COOP, sorted_final_values, 'varnames', {'FinalValue'}, 'display', 'off');

% Tukey HSD (multiple comparison)
[c, m, h, gnames] = multcompare(stats, 'CType', 'hsd', 'Display', 'off');


% ---- Plot graph ----
% Combine data for 'No Obstruction' and 'With Obstruction', calculate mean and standard deviation
unique_final_values = categories(FinalValue);
unique_final_values_num = str2double(unique_final_values);

COOP_mean = zeros(length(unique_final_values),1);
COOP_sd = zeros(length(unique_final_values),1);

for i = 1:length(unique_final_values)
    fv = unique_final_values{i};
    idx = FinalValue == fv;
    COOP_mean(i) = mean(sorted_COOP(idx));
    COOP_sd(i) = std(sorted_COOP(idx));
end

% Create bar plot
figure;
bar(unique_final_values_num, COOP_mean, 'FaceColor', 'b', 'BarWidth', 0.5);
hold on;
errorbar(unique_final_values_num, COOP_mean, COOP_sd, '.', 'Color', 'k');
xlabel('Final Value (Stretch Condition)');
ylabel('Mean COOP');
ylim([0, 1]);
title('COOP for Different Final Values');
grid on;

% Display significant differences (only for p < 0.05)
% Set starting position and increment for y-axis for significance lines
y_base = max(COOP_mean + COOP_sd) + 0.1; % Set the starting position of the line
y_increment = 0.05; % Increment for the height of the line

% Initialize counter to record the current number of lines
line_counter = 0;

for i = 1:length(unique_final_values)
    for j = i+1:length(unique_final_values)
        % Get group numbers
        group1 = find(strcmp(gnames, ['FinalValue=' unique_final_values{i}]));
        group2 = find(strcmp(gnames, ['FinalValue=' unique_final_values{j}]));
        
        % Find pairwise comparisons
        comp_idx = ( (c(:,1) == group1 & c(:,2) == group2) | (c(:,1) == group2 & c(:,2) == group1) );
        if any(comp_idx)
            p_value = c(comp_idx, 6);
            if p_value < 0.05
                % Plot significant differences
                x1 = unique_final_values_num(i) + 0.1;
                x2 = unique_final_values_num(j) - 0.1;
                y = y_base + line_counter * y_increment;
                plot([x1, x2], [y, y], '-k', 'LineWidth', 1.5);
                text(mean([x1, x2]), y + 0.02, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
                % Increment counter
                line_counter = line_counter + 1;
            end
        end
    end
end


hold off;

% ---- Display results ----
disp('ANOVA results for COOP:');
disp(tbl);
disp('Tukey HSD multiple comparison results:');
disp(array2table(c, 'VariableNames', {'Group1','Group2','LowerLimit','MeanDiff','UpperLimit','pValue'}));


%% Function definitions

% Function: Calculate OOP and fiber orientation grid
function [OOP, grid] = calculate_OOP_and_fiber_orientation_grid(filename)
    load(filename, 'net_res');
    grid_size = 10; % Number of grids
    grid = zeros(grid_size, grid_size, 2); % Store average direction vectors
    total_vectors = []; % All fiber direction vectors

    for sim_idx = 1:length(net_res)
        net = net_res{1, sim_idx}{end}; % Final step
        % Collect vectors for each grid
        [grid_vectors, vectors] = collect_grid_vectors(net, grid_size);
        grid = grid + grid_vectors;
        total_vectors = [total_vectors, vectors];
    end

    % Calculate OOP
    OOP = calculate_OOP(total_vectors);
end

% Function: Calculate COOP
function COOP = calculate_COOP(OOP1, OOP2, grid1, grid2)
    N = size(grid1, 1) * size(grid1, 2); % Number of grids
    T_sum = zeros(2, 2);

    for xi = 1:size(grid1, 1)
        for yi = 1:size(grid1, 2)
            % Get grid vectors
            v1 = squeeze(grid1(xi, yi, :));
            v2 = squeeze(grid2(xi, yi, :));

            if any(v1) && any(v2) % Calculate only non-zero vectors
                % Extend vectors to 3D (set z component to 0)
                v1_3d = [v1; 0];
                v2_3d = [v2; 0];

                % Calculate cos(θ) and sin(θ)
                f_x = dot(v1, v2); % cos(θ)
                f_y = norm(cross(v1_3d, v2_3d)); % sin(θ)

                % Calculate f_i
                f_i = [f_x; f_y];

                % Calculate T_i
                T_i = 2 * (f_i * f_i') - eye(2);
                T_sum = T_sum + T_i;
            end
        end
    end

    % Calculate average tensor
    T_mean = T_sum / N;
    COOP = max(eig(T_mean)); % Maximum eigenvalue
end

% Function: Collect fiber orientation vectors for each grid
function [grid_vectors, total_vectors] = collect_grid_vectors(net, grid_size)
    % Initialize grid edges
    all_x = [];
    all_y = [];
    for j = 1:length(net)
        net_temp = net{j};
        if isempty(net_temp)
            continue;
        end
        % Consider the case where fiber data is a 3D array
        if size(net_temp, 1) >= 2
            x_coords = net_temp(1, :, :);
            y_coords = net_temp(2, :, :);
            % Convert to 2D array (reshape multidimensional array to vector)
            x_coords = x_coords(:);
            y_coords = y_coords(:);
            all_x = [all_x; x_coords];
            all_y = [all_y; y_coords];
        else
            warning('Unexpected net_temp dimensions in net{%d}', j);
        end
    end

    % Define grid edges
    x_min = min(all_x); x_max = max(all_x);
    y_min = min(all_y); y_max = max(all_y);
    x_edges = linspace(x_min, x_max, grid_size+1);
    y_edges = linspace(y_min, y_max, grid_size+1);

    % Collect vectors for each grid
    grid_vectors = zeros(grid_size, grid_size, 2); % Store average vectors
    total_vectors = [];

    for xi = 1:grid_size
        for yi = 1:grid_size
            % Range of grid cell
            x_start = x_edges(xi); x_end = x_edges(xi+1);
            y_start = y_edges(yi); y_end = y_edges(yi+1);

            % Vectors within the grid cell
            cell_vectors = [];

            for j = 1:length(net)
                net_temp = net{j};
                if isempty(net_temp)
                    continue;
                end

                for idx = 1:size(net_temp, 2)-1
                    % Start and end points of the fiber
                    if size(net_temp, 1) >= 2
                        p1 = net_temp(:, idx);
                        p2 = net_temp(:, idx+1);
                        mid_point = (p1 + p2) / 2;

                        % Determine if the midpoint is within the grid
                        if mid_point(1) >= x_start && mid_point(1) < x_end && ...
                           mid_point(2) >= y_start && mid_point(2) < y_end
                            % Calculate vector
                            vec = p2 - p1;
                            vec = vec / norm(vec); % Normalize
                            cell_vectors = [cell_vectors, vec];
                        end
                    else
                        warning('Unexpected net_temp dimensions in net{%d}', j);
                    end
                end
            end

            % Calculate average vector for the grid
            if ~isempty(cell_vectors)
                avg_vec = mean(cell_vectors, 2);
                avg_vec = avg_vec / norm(avg_vec);
                grid_vectors(xi, yi, :) = avg_vec;
                total_vectors = [total_vectors, avg_vec];
            end
        end
    end
end

% Function: Calculate OOP
function OOP = calculate_OOP(vectors)
    N = size(vectors, 2);
    if N > 0
        T_sum = zeros(2, 2);
        for i = 1:N
            v = vectors(:, i);
            T_i = 2 * (v * v') - eye(2);
            T_sum = T_sum + T_i;
        end
        T_mean = T_sum / N;
        OOP = max(eig(T_mean));
    else
        OOP = NaN;
    end
end