function max_force = calculate_last_traction_force(F_store)
    % Function to calculate the maximum traction force
    % F_store: Force data (Nx2 matrix, [F_x, F_y])
    % max_force: Maximum traction force

    % Calculate the magnitude of the traction force at each time point
    forces = sqrt(F_store(:, :, 1).^2 + F_store(:, :, 2).^2); % Resultant force
    % max_force = max(forces(:)); % Get the maximum value
    % Use last time step force
    % Get the maximum absolute value
    max_force = max(forces(end, :));
end
