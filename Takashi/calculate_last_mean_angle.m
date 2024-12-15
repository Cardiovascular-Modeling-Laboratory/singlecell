function angle_deg = calculate_last_mean_angle_v3(net_res)
    % Function to calculate OOP (returns the final mean angle as a scalar)
    grid_size = 10; % Number of grids
    grid = zeros(grid_size, grid_size, 2); 
    total_vectors = []; % All fiber direction vectors

    for sim_idx = 1:length(net_res)
        net = net_res{1, sim_idx}{end}; % Final step
        [grid_vectors, vectors] = collect_grid_vectors(net, grid_size);
        grid = grid + grid_vectors;
        total_vectors = [total_vectors, vectors];
    end

    OOP = calculate_angle_from_vectors(total_vectors);
    angle_deg = rad2deg(OOP);
end

function [grid_vectors, total_vectors] = collect_grid_vectors(net, grid_size)
    all_x = [];
    all_y = [];
    for j = 1:length(net)
        net_temp = net{j};
        if isempty(net_temp)
            continue;
        end
        if size(net_temp, 1) >= 2
            x_coords = net_temp(1, :, :);
            y_coords = net_temp(2, :, :);
            x_coords = x_coords(:);
            y_coords = y_coords(:);
            all_x = [all_x; x_coords];
            all_y = [all_y; y_coords];
        else
            warning('Unexpected net_temp dimensions in net{%d}', j);
        end
    end

    x_min = min(all_x); x_max = max(all_x);
    y_min = min(all_y); y_max = max(all_y);
    x_edges = linspace(x_min, x_max, grid_size+1);
    y_edges = linspace(y_min, y_max, grid_size+1);

    grid_vectors = zeros(grid_size, grid_size, 2);
    total_vectors = [];
    
    for xi = 1:grid_size
        for yi = 1:grid_size
            x_start = x_edges(xi); x_end = x_edges(xi+1);
            y_start = y_edges(yi); y_end = y_edges(yi+1);
            
            cell_vectors = [];
            
            for j = 1:length(net)
                net_temp = net{j};
                if isempty(net_temp)
                    continue;
                end
                
                for idx = 1:size(net_temp, 2)-1
                    if size(net_temp, 1) >= 2
                        p1 = net_temp(:, idx);
                        p2 = net_temp(:, idx+1);
                        mid_point = (p1 + p2) / 2;
                        if mid_point(1) >= x_start && mid_point(1) < x_end && ...
                           mid_point(2) >= y_start && mid_point(2) < y_end
                            vec = p2 - p1;
                            vec = vec / norm(vec); % Normalize
                            
                            % Map angle to 0-90 degrees
                            angle = atan2(vec(2), vec(1));
                            angle_deg = rad2deg(angle); % [-180,180)
                            if angle_deg < 0
                                angle_deg = angle_deg + 180; % [0,180)
                            end
                            if angle_deg > 90
                                angle_deg = 180 - angle_deg; % Fit within [0,90]
                            end
                            
                            vec_new = [cosd(angle_deg); sind(angle_deg)];
                            cell_vectors = [cell_vectors, vec_new];
                        end
                    else
                        warning('Unexpected net_temp dimensions in net{%d}', j);
                    end
                end
            end
            
            if ~isempty(cell_vectors)
                avg_vec = mean(cell_vectors, 2);
                avg_vec = avg_vec / norm(avg_vec);

                % Optionally unify avg_vec to 0-90 degrees (can be omitted if unnecessary)
                angle_avg = atan2(avg_vec(2), avg_vec(1));
                angle_avg_deg = rad2deg(angle_avg);
                if angle_avg_deg < 0
                    angle_avg_deg = angle_avg_deg + 180;
                end
                if angle_avg_deg > 90
                    angle_avg_deg = 180 - angle_avg_deg;
                end
                avg_vec = [cosd(angle_avg_deg); sind(angle_avg_deg)];
                
                grid_vectors(xi, yi, :) = avg_vec;
                total_vectors = [total_vectors, avg_vec];
            end
        end
    end
end

function angle = calculate_angle_from_vectors(vectors)
    N = size(vectors, 2);
    if N > 0
        T_sum = zeros(2, 2);
        for i = 1:N
            v = vectors(:, i);
            T_i = 2 * (v * v') - eye(2);
            T_sum = T_sum + T_i;
        end
        T_mean = T_sum / N;
        
        [V, D] = eig(T_mean);
        [~, max_idx] = max(diag(D));
        principal_vector = V(:, max_idx);
        angle = atan2(principal_vector(2), principal_vector(1));
        
        % Fit principal_vector within 0-90 degrees
        angle_deg = rad2deg(angle);
        if angle_deg < 0
            angle_deg = angle_deg + 180;
        end
        if angle_deg > 90
            angle_deg = 180 - angle_deg;
        end
        angle = deg2rad(angle_deg);
        
        %% Debug plot %%
        figure; 
        % 1. Visualize distribution of all vectors
        subplot(1,3,1);
        quiver(zeros(1,N), zeros(1,N), vectors(1,:), vectors(2,:), 0, 'b');
        hold on;
        % Plot principal direction vector in red
        quiver(0,0,cosd(angle_deg), sind(angle_deg),0,'r','LineWidth',2);
        axis equal; grid on;
        title('All fiber direction vectors (0-90 deg)');
        xlabel('X'); ylabel('Y');
        
        % 2. Polar histogram of angle distribution
        subplot(1,3,2);
        all_angles = atan2(vectors(2,:), vectors(1,:));
        all_angles_deg = rad2deg(all_angles);
        % Normalize to 0-90 degrees
        all_angles_deg(all_angles_deg<0) = all_angles_deg(all_angles_deg<0) + 180;
        mask = all_angles_deg > 90;
        all_angles_deg(mask) = 180 - all_angles_deg(mask);
        polarhistogram(deg2rad(all_angles_deg), 18);
        title('Angle distribution of fibers (0-90 deg)');
        
        % 3. Visualize T_mean tensor ellipse
        subplot(1,3,3);
        [V_plot, D_plot] = eig(T_mean);
        th = linspace(0,2*pi,100);
        unit_circle = [cos(th); sin(th)];
        scale_matrix = sqrt(D_plot);
        ellipse_points = V_plot * scale_matrix * unit_circle;
        plot(ellipse_points(1,:), ellipse_points(2,:),'k-','LineWidth',1.5);
        hold on;
        quiver(0,0,cosd(angle_deg), sind(angle_deg), 0,'r','LineWidth',2);
        axis equal; grid on;
        title('T mean ellipse and principal direction');
        xlabel('X'); ylabel('Y');
        
    else
        angle = NaN;
    end
end
