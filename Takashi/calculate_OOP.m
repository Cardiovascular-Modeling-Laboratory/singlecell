function OOP = calculate_OOP(net_res)
    % Function to calculate OOP
    % net_res: Fiber data from simulation results
    % OOP: Orientation Order Parameter (scalar value)

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
    OOP = calculate_OOP_from_vectors(total_vectors);
end

%% Function: Collect fiber direction vectors for each grid
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

%% Function: Calculate OOP
function OOP = calculate_OOP_from_vectors(vectors)
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
