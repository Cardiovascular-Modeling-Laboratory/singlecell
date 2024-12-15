% Script to plot necessary data

dataDir = pwd; % Current working directory
% Get .mat files in the directory (change file pattern to include results with n=3)
filePattern = fullfile(dataDir, 'Takashi', 'results', 'fig3_v2_n=*_stretch_*.mat');
resultFiles = dir(filePattern);

resultFiles

% Variables for saving data (changed to accumulate data)
stretch_factors = [];
mean_angles_all = [];
OOP_values_all = [];
traction_forces_all = [];

% Collect data
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);
    
    % Extract stretch_factor and n from the filename
    tokens = regexp(resultFiles(k).name, 'fig3_v2_n=(\d+)_stretch_(\d+)', 'tokens');
    if ~isempty(tokens)
        n_sim = str2double(tokens{1}{1});
        stretch_factor = str2double(tokens{1}{2});
    else
        disp(['Could not parse n or stretch_factor from filename: ', resultFiles(k).name]);
        continue; % Move to the next file
    end
    
    % Load data
    load(filename, 'net_res', 'F_res', 'F_store');  % Load necessary data

    % ---- Calculate and save necessary data ----
    % Calculate mean angles and OOP from multiple simulation results
    mean_angles = calculate_last_mean_angle(net_res);      % Size: [n_sim x 1]
    OOP_values = calculate_OOP(net_res);        % Size: [n_sim x 1]
    traction_forces = calculate_last_traction_force(F_store); % Size: [n_sim x 1]
    
    % Save results
    stretch_factors = [stretch_factors; stretch_factor];
    mean_angles_all = [mean_angles_all; mean_angles];
    OOP_values_all = [OOP_values_all; OOP_values];
    traction_forces_all = [traction_forces_all; traction_forces];
end

% Get unique stretch factors
unique_stretch_factors = unique(stretch_factors);

% ---- Add n number limit to data collection ----
filtered_stretch_factors = [];
filtered_mean_angles_all = [];
filtered_OOP_values_all = [];
filtered_traction_forces_all = [];

for i = 1:length(unique_stretch_factors)
    sf = unique_stretch_factors(i);
    idx = find(stretch_factors == sf); % Get indices that meet the condition
    
    % Limit n number to 3 (select the first 3, random sampling can be used if needed)
    if length(idx) > 3
        idx = idx(1:3); % Select the first 3
        % idx = idx(randperm(length(idx), 3)); % Use this for random selection of 3
    end
    
    % Save data after limiting
    filtered_stretch_factors = [filtered_stretch_factors; stretch_factors(idx)];
    filtered_mean_angles_all = [filtered_mean_angles_all; mean_angles_all(idx)];
    filtered_OOP_values_all = [filtered_OOP_values_all; OOP_values_all(idx)];
    filtered_traction_forces_all = [filtered_traction_forces_all; traction_forces_all(idx)];
    
    % Display n number (always 3 after limiting)
    n_count = length(idx); % n number after limiting
    disp(['Stretch factor: ', num2str(sf), ', n: ', num2str(n_count)]);
end

% ---- Calculate mean and standard deviation (after limiting) ----
mean_angles_mean = zeros(length(unique_stretch_factors),1);
mean_angles_sd = zeros(length(unique_stretch_factors),1);
OOP_values_mean = zeros(length(unique_stretch_factors),1);
OOP_values_sd = zeros(length(unique_stretch_factors),1);
traction_forces_mean = zeros(length(unique_stretch_factors),1);
traction_forces_sd = zeros(length(unique_stretch_factors),1);

for i = 1:length(unique_stretch_factors)
    sf = unique_stretch_factors(i);
    idx = filtered_stretch_factors == sf; % Apply condition to filtered data
    
    mean_angles_mean(i) = mean(filtered_mean_angles_all(idx));
    mean_angles_sd(i) = std(filtered_mean_angles_all(idx));
    OOP_values_mean(i) = mean(filtered_OOP_values_all(idx));
    OOP_values_sd(i) = std(filtered_OOP_values_all(idx));
    traction_forces_mean(i) = mean(filtered_traction_forces_all(idx));
    traction_forces_sd(i) = std(filtered_traction_forces_all(idx));
end

% ---- Display n number in data collection ----
for i = 1:length(unique_stretch_factors)
    sf = unique_stretch_factors(i);
    idx = stretch_factors == sf;
    n_count = sum(idx); % Number of data that meet the condition
    disp(['Stretch factor: ', num2str(sf), ', n: ', num2str(n_count)]);
end

% ANOVA and post-hoc test for Mean Angle
[p_angle, tbl_angle, stats_angle] = anova1(filtered_mean_angles_all, filtered_stretch_factors, 'off');
[c_angle, m_angle, h_angle, gnames_angle] = multcompare(stats_angle, 'Display', 'off');

% ANOVA and post-hoc test for OOP
[p_OOP, tbl_OOP, stats_OOP] = anova1(filtered_OOP_values_all, filtered_stretch_factors, 'off');
[c_OOP, m_OOP, h_OOP, gnames_OOP] = multcompare(stats_OOP, 'Display', 'off');

% ANOVA and post-hoc test for Traction Force
[p_traction, tbl_traction, stats_traction] = anova1(filtered_traction_forces_all, filtered_stretch_factors, 'off');
[c_traction, m_traction, h_traction, gnames_traction] = multcompare(stats_traction, 'Display', 'off');

debug_show_asterisks = false;
% ---- Plot graphs ----
figure;

% Graph 1: Mean Angle vs Stretch Factor
subplot(1, 3, 1);
bar(unique_stretch_factors, mean_angles_mean);
hold on;
errorbar(unique_stretch_factors, mean_angles_mean, mean_angles_sd, 'LineStyle', 'none', 'LineWidth', 2); % Change thickness
xlabel('Stretch Factor');
ylabel('Mean Angle (degrees)');
ylim([0, 360]); % Set angle range to 0-180

sig_positions = max(mean_angles_mean + mean_angles_sd) + 5;
significant_pairs = c_angle(c_angle(:,6) < 0.05, 1:2); % Pairs with p<0.05
offset = 0;
for i = 1:size(significant_pairs, 1)
    group1 = significant_pairs(i,1);
    group2 = significant_pairs(i,2);
    x1 = unique_stretch_factors(group1);
    x2 = unique_stretch_factors(group2);
    y = sig_positions + offset;
    plot([x1, x1, x2, x2], [y, y+1, y+1, y], '-k', 'LineWidth', 1);
    text(mean([x1, x2]), y+1.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    offset = offset + 5; % Adjust height to avoid overlapping lines
end

% Graph 2: OOP vs Stretch Factor
subplot(1, 3, 2);
bar(unique_stretch_factors, OOP_values_mean);
hold on;
errorbar(unique_stretch_factors, OOP_values_mean, OOP_values_sd, 'LineStyle', 'none', 'LineWidth', 2); % Change thickness
xlabel('Stretch Factor');
ylabel('OOP');
title('OOP');
ylim([0, 1]); % Set OOP range to 0-1

sig_positions = max(OOP_values_mean + OOP_values_sd) + 0.05;
offset = 0;
significant_pairs = c_OOP(c_OOP(:,6) < 0.05, 1:2);
for i = 1:size(significant_pairs, 1)
    group1 = significant_pairs(i,1);
    group2 = significant_pairs(i,2);
    x1 = unique_stretch_factors(group1);
    x2 = unique_stretch_factors(group2);
    y = sig_positions + offset;
    plot([x1, x1, x2, x2], [y, y+0.02, y+0.02, y], '-k', 'LineWidth', 1);
    text(mean([x1, x2]), y+0.03, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    offset = offset + 0.05; % Adjust height to avoid overlapping lines
end

% Graph 3: Traction Force vs Stretch Factor
subplot(1, 3, 3);
bar(unique_stretch_factors, traction_forces_mean);
hold on;
errorbar(unique_stretch_factors, traction_forces_mean, traction_forces_sd, 'LineStyle', 'none', 'LineWidth', 2); % Change thickness
xlabel('Stretch Factor');
ylabel('Traction Force (N)');

sig_positions = max(traction_forces_mean + traction_forces_sd) + 0.05;
offset = 0;
significant_pairs = c_traction(c_traction(:,6) < 0.05, 1:2);
for i = 1:size(significant_pairs, 1)
    group1 = significant_pairs(i,1);
    group2 = significant_pairs(i,2);
    x1 = unique_stretch_factors(group1);
    x2 = unique_stretch_factors(group2);
    y = sig_positions + offset;
    plot([x1, x1, x2, x2], [y, y+0.02, y+0.02, y], '-k', 'LineWidth', 1);
    text(mean([x1, x2]), y+0.03, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    offset = offset + 0.05; % Adjust height to avoid overlapping lines
end

% Adjust overall layout
sgtitle('Stretch Results Analysis');

% ---- Statistical tests by ANOVA ----
% Perform one-way ANOVA for each variable
[p_angle, tbl_angle] = anova1(mean_angles_all, stretch_factors, 'off');
[p_OOP, tbl_OOP] = anova1(OOP_values_all, stretch_factors, 'off');
[p_traction, tbl_traction] = anova1(traction_forces_all, stretch_factors, 'off');

% Display results
disp('ANOVA for Mean Angles:');
disp(tbl_angle);
disp(['p-value: ', num2str(p_angle)]);

disp('ANOVA for OOP:');
disp(tbl_OOP);
disp(['p-value: ', num2str(p_OOP)]);

disp('ANOVA for Traction Forces:');
disp(tbl_traction);
disp(['p-value: ', num2str(p_traction)]);



