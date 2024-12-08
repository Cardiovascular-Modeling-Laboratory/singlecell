% filepath: /Users/inagakit/Documents/UCIrvine/AnnaGrosberg/singlecell/Takashi/calculate_last_mean_angle_v3.m
function mean_angles_grid = calculate_last_mean_angle_v3(net_res)
% Calculate mean angles per grid cell using the mean tensor from OOP
% net_res: Simulation result containing fiber data
% grid_size: Number of grid divisions along each axis
% mean_angles_grid: Grid of mean angles in degrees

% Initialize grid for mean angles
grid_size = 10; % グリッドの数
mean_angles_grid = NaN(grid_size, grid_size);

% Use the data from the last step of net_res{1,1}
net1 = net_res{1, 1};
net = net1{end};

% Collect all coordinates to define grid edges
all_x = [];
all_y = [];
for i = 1:length(net)
    net_temp = net{i};
    if isempty(net_temp)
        continue;
    end
    % Extract coordinates
    x_coords = net_temp(1, :, :);
    y_coords = net_temp(2, :, :);
    x_coords = x_coords(:);
    y_coords = y_coords(:);
    all_x = [all_x; x_coords];
    all_y = [all_y; y_coords];
end

% Define grid edges
x_min = min(all_x); x_max = max(all_x);
y_min = min(all_y); y_max = max(all_y);
x_edges = linspace(x_min, x_max, grid_size+1);
y_edges = linspace(y_min, y_max, grid_size+1);

% Initialize cell arrays to store direction vectors per grid cell
grid_vectors = cell(grid_size, grid_size);

% Loop through fibers and assign direction vectors to grid cells
for i = 1:length(net)
    net_temp = net{i};
    if isempty(net_temp)
        continue;
    end
    for j = 1:size(net_temp, 3)
        fiber = net_temp(:, :, j);
        for idx = 1:size(fiber, 2) - 1
            p1 = fiber(:, idx);
            p2 = fiber(:, idx + 1);
            v = p2 - p1;
            v = v / norm(v); % Normalize

            % Determine grid cell for p1
            xi = find(x_edges <= p1(1), 1, 'last');
            yi = find(y_edges <= p1(2), 1, 'last');
            xi = min(xi, grid_size);
            yi = min(yi, grid_size);

            % Store vector in corresponding grid cell
            grid_vectors{xi, yi}{end+1} = v;
        end
    end
end

% Compute mean tensor and mean angle per grid cell
for xi = 1:grid_size
    for yi = 1:grid_size
        vectors = grid_vectors{xi, yi};
        if ~isempty(vectors)
            T_sum = zeros(2, 2);
            for k = 1:length(vectors)
                v = vectors{k};
                T_i = 2 * (v * v') - eye(2);
                T_sum = T_sum + T_i;
            end
            T_mean = T_sum / length(vectors);
            [vecs, vals] = eig(T_mean);
            [~, idx] = max(diag(vals));
            principal_direction = vecs(:, idx);
            mean_angle = atan2d(principal_direction(2), principal_direction(1));
            mean_angles_grid(xi, yi) = mean_angle;
        end
    end
end
end