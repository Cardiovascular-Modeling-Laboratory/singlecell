% filepath: /Users/inagakit/Documents/UCIrvine/AnnaGrosberg/singlecell/Takashi/calculate_last_mean_angle_v2.m
function mean_angle = calculate_last_mean_angle_v2(net_res)
% Mean angle calculation using the mean tensor from OOP
% net_res: Simulation result containing fiber data
% mean_angle: Mean angle in degrees

% Initialize tensor sum
T_sum = zeros(2, 2);
N = 0;

% Use the data from the last step of net_res{1,1}
net1 = net_res{1, 1};
net = net1{end};

% Loop through fibers
for i = 1:length(net)
    net_temp = net{i};
    if isempty(net_temp)
        continue;
    end

    % Process each fiber in net_temp
    for j = 1:size(net_temp, 3)
        fiber = net_temp(:, :, j);

        % Calculate direction vectors along the fiber
        for idx = 1:size(fiber, 2) - 1
            p1 = fiber(:, idx);
            p2 = fiber(:, idx + 1);
            v = p2 - p1;
            v = v / norm(v); % Normalize
            T_i = 2 * (v * v') - eye(2);
            T_sum = T_sum + T_i;
            N = N + 1;
        end
    end
end

% Calculate mean tensor
T_mean = T_sum / N;

% Compute eigenvalues and eigenvectors
[vecs, vals] = eig(T_mean);

% Find the principal direction (eigenvector with the maximum eigenvalue)
[max_eigval, idx] = max(diag(vals));
principal_direction = vecs(:, idx);

% Compute mean angle from the principal direction
mean_angle = atan2d(principal_direction(2), principal_direction(1));

% ...existing code...
end