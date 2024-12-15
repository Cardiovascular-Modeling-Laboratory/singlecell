% Script to plot necessary data

% Get .mat files in the directory (change file pattern to include results with n=3)
dataDir = pwd;
filePattern = fullfile(dataDir, 'Takashi', 'results', 'fig9_v2_n=*_final*.mat');
resultFiles = dir(filePattern);

% Variables for saving data (changed to accumulate data)
final_values = [];
mean_angles_all = [];
OOP_values_all = [];
traction_forces_all = [];
nucleus_obstruction_all = [];

% Collect data
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);

    % Extract n, final_value, and nucleus_obstruction from the filename
    tokens = regexp(resultFiles(k).name, 'fig9_v2_n=(\d+)_final(\d+)_nuc_rel_vec_(\d+)', 'tokens');
    if ~isempty(tokens)
        n_sim = str2double(tokens{1}{1});
        final_value = str2double(tokens{1}{2});
        nuc_obstruction = str2double(tokens{1}{3});
    else
        disp(['Could not parse n, final value, or nucleus obstruction from filename: ', resultFiles(k).name]);
        continue; % Move to the next file
    end

    % Load data
    load(filename, 'net_res', 'F_res', 'F_store');  % Load necessary data

    % ---- Calculate and save necessary data ----
    % Calculate mean angle and OOP from multiple simulation results
    mean_angles = calculate_last_mean_angle(net_res);      % Size: [n_sim x 1]
    OOP_values = calculate_OOP(net_res);                   % Size: [n_sim x 1]
    traction_forces = calculate_last_traction_force(F_store); % Size: [n_sim x 1]

    % Save results
    final_values = [final_values; final_value];
    mean_angles_all = [mean_angles_all; mean_angles];
    OOP_values_all = [OOP_values_all; OOP_values];
    traction_forces_all = [traction_forces_all; traction_forces];
    nucleus_obstruction_all = [nucleus_obstruction_all; nuc_obstruction];
end

% Sort data
[sorted_final_values, sort_idx] = sort(final_values);
sorted_mean_angles = mean_angles_all(sort_idx);
sorted_OOP_values = OOP_values_all(sort_idx);
sorted_traction_forces = traction_forces_all(sort_idx);
sorted_nucleus_obstruction = nucleus_obstruction_all(sort_idx);

% Split data by presence or absence of nucleus obstruction
no_obstruction_idx = sorted_nucleus_obstruction == 0;
with_obstruction_idx = sorted_nucleus_obstruction == 1;

% Get unique final values
unique_final_values = unique(sorted_final_values);

% Split data by each final value and calculate mean and standard deviation for presence or absence of nucleus
mean_angles_mean_no = zeros(length(unique_final_values),1);
mean_angles_sd_no = zeros(length(unique_final_values),1);
mean_angles_mean_with = zeros(length(unique_final_values),1);
mean_angles_sd_with = zeros(length(unique_final_values),1);

OOP_values_mean_no = zeros(length(unique_final_values),1);
OOP_values_sd_no = zeros(length(unique_final_values),1);
OOP_values_mean_with = zeros(length(unique_final_values),1);
OOP_values_sd_with = zeros(length(unique_final_values),1);

traction_forces_mean_no = zeros(length(unique_final_values),1);
traction_forces_sd_no = zeros(length(unique_final_values),1);
traction_forces_mean_with = zeros(length(unique_final_values),1);
traction_forces_sd_with = zeros(length(unique_final_values),1);

for i = 1:length(unique_final_values)
    fv = unique_final_values(i);
    % No Obstruction
    idx_no = sorted_final_values == fv & no_obstruction_idx;
    mean_angles_mean_no(i) = mean(sorted_mean_angles(idx_no));
    mean_angles_sd_no(i) = std(sorted_mean_angles(idx_no));
    OOP_values_mean_no(i) = mean(sorted_OOP_values(idx_no));
    OOP_values_sd_no(i) = std(sorted_OOP_values(idx_no));
    traction_forces_mean_no(i) = mean(sorted_traction_forces(idx_no));
    traction_forces_sd_no(i) = std(sorted_traction_forces(idx_no));
    % With Obstruction
    idx_with = sorted_final_values == fv & with_obstruction_idx;
    mean_angles_mean_with(i) = mean(sorted_mean_angles(idx_with));
    mean_angles_sd_with(i) = std(sorted_mean_angles(idx_with));
    OOP_values_mean_with(i) = mean(sorted_OOP_values(idx_with));
    OOP_values_sd_with(i) = std(sorted_OOP_values(idx_with));
    traction_forces_mean_with(i) = mean(sorted_traction_forces(idx_with));
    traction_forces_sd_with(i) = std(sorted_traction_forces(idx_with));
end

% Adjust x-axis positions
x_no_obstruction = unique_final_values;
x_with_obstruction = unique_final_values + 0.35;

% Debug flag
debug_show_asterisks = true;

% ---- Plot graphs ----
figure;

% Graph 1: Mean Angle vs Final Value
subplot(1, 3, 1);
hold on;
bar1 = bar(x_no_obstruction, mean_angles_mean_no, 'FaceColor', 'b', 'BarWidth', 0.3);
errorbar(x_no_obstruction, mean_angles_mean_no, mean_angles_sd_no, '.', 'Color', 'k');
bar2 = bar(x_with_obstruction, mean_angles_mean_with, 'FaceColor', 'r', 'BarWidth', 0.3);
errorbar(x_with_obstruction, mean_angles_mean_with, mean_angles_sd_with, '.', 'Color', 'k');
xlabel('Final Value');
ylabel('Mean Angle (degrees)');
ylim([0, 90]); % Set angle range to 0-90
legend([bar1, bar2], {'No Obstruction', 'With Obstruction'});
title('Mean Angle');

% Display significance based on ANOVA results
sig_positions = max([mean_angles_mean_no + mean_angles_sd_no, mean_angles_mean_with + mean_angles_sd_with], [], 2) + 5;
for i = 1:length(unique_final_values)
    idx_no = sorted_final_values == unique_final_values(i) & sorted_nucleus_obstruction == 0;
    idx_with = sorted_final_values == unique_final_values(i) & sorted_nucleus_obstruction == 1;
    [~, p] = ttest2(sorted_mean_angles(idx_no), sorted_mean_angles(idx_with));
    if p < 0.05 || debug_show_asterisks
        plot([x_no_obstruction(i), x_with_obstruction(i)], [sig_positions(i), sig_positions(i)], '-k', 'LineWidth', 1.5);
        text(mean([x_no_obstruction(i), x_with_obstruction(i)]), sig_positions(i) + 2, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end

% Graph 2: OOP vs Final Value
subplot(1, 3, 2);
hold on;
bar1 = bar(x_no_obstruction, OOP_values_mean_no, 'FaceColor', 'b', 'BarWidth', 0.3);
errorbar(x_no_obstruction, OOP_values_mean_no, OOP_values_sd_no, '.', 'Color', 'k');
bar2 = bar(x_with_obstruction, OOP_values_mean_with, 'FaceColor', 'r', 'BarWidth', 0.3);
errorbar(x_with_obstruction, OOP_values_mean_with, OOP_values_sd_with, '.', 'Color', 'k');
xlabel('Final Value');
ylabel('OOP');
ylim([0, 1]); % Set OOP range to 0-1
legend([bar1, bar2], {'No Obstruction', 'With Obstruction'});
title('OOP');

% Display significance based on ANOVA results
sig_positions = max([OOP_values_mean_no + OOP_values_sd_no, OOP_values_mean_with + OOP_values_sd_with], [], 2) + 0.05;
for i = 1:length(unique_final_values)
    idx_no = sorted_final_values == unique_final_values(i) & sorted_nucleus_obstruction == 0;
    idx_with = sorted_final_values == unique_final_values(i) & sorted_nucleus_obstruction == 1;
    [~, p] = ttest2(sorted_OOP_values(idx_no), sorted_OOP_values(idx_with));
    if p < 0.05 || debug_show_asterisks
        plot([x_no_obstruction(i), x_with_obstruction(i)], [sig_positions(i), sig_positions(i)], '-k', 'LineWidth', 1.5);
        text(mean([x_no_obstruction(i), x_with_obstruction(i)]), sig_positions(i) + 0.02, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end

% Graph 3: Traction Force vs Final Value
subplot(1, 3, 3);
hold on;
bar1 = bar(x_no_obstruction, traction_forces_mean_no, 'FaceColor', 'b', 'BarWidth', 0.3);
errorbar(x_no_obstruction, traction_forces_mean_no, traction_forces_sd_no, '.', 'Color', 'k');
bar2 = bar(x_with_obstruction, traction_forces_mean_with, 'FaceColor', 'r', 'BarWidth', 0.3);
errorbar(x_with_obstruction, traction_forces_mean_with, traction_forces_sd_with, '.', 'Color', 'k');
xlabel('Final Value');
ylabel('Traction Force (N)');
legend([bar1, bar2], {'No Obstruction', 'With Obstruction'});
title('Traction Force');

% Display significance based on ANOVA results
sig_positions = max([traction_forces_mean_no + traction_forces_sd_no, traction_forces_mean_with + traction_forces_sd_with], [], 2) + 0.05;
for i = 1:length(unique_final_values)
    idx_no = sorted_final_values == unique_final_values(i) & sorted_nucleus_obstruction == 0;
    idx_with = sorted_final_values == unique_final_values(i) & sorted_nucleus_obstruction == 1;
    [~, p] = ttest2(sorted_traction_forces(idx_no), sorted_traction_forces(idx_with));
    if p < 0.05 || debug_show_asterisks
        plot([x_no_obstruction(i), x_with_obstruction(i)], [sig_positions(i), sig_positions(i)], '-k', 'LineWidth', 1.5);
        text(mean([x_no_obstruction(i), x_with_obstruction(i)]), sig_positions(i) + 0.02, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end

% Adjust overall layout
sgtitle('Stretch Results Analysis with and without Nucleus Obstruction');

% ---- Statistical test by ANOVA ----
% Perform two-way ANOVA
group = {sorted_final_values, sorted_nucleus_obstruction};

% ANOVA for Mean Angle
[p_angle, tbl_angle, stats_angle] = anovan(sorted_mean_angles, group, 'model', 'interaction', 'varnames', {'FinalValue', 'NucleusObstruction'});

% ANOVA for OOP
[p_OOP, tbl_OOP, stats_OOP] = anovan(sorted_OOP_values, group, 'model', 'interaction', 'varnames', {'FinalValue', 'NucleusObstruction'});

% ANOVA for Traction Force
[p_traction, tbl_traction, stats_traction] = anovan(sorted_traction_forces, group, 'model', 'interaction', 'varnames', {'FinalValue', 'NucleusObstruction'});

% Display results
disp('ANOVA for Mean Angles:');
disp(tbl_angle);
disp(['p-values: ', num2str(p_angle')]);

disp('ANOVA for OOP:');
disp(tbl_OOP);
disp(['p-values: ', num2str(p_OOP')]);

disp('ANOVA for Traction Forces:');
disp(tbl_traction);
disp(['p-values: ', num2str(p_traction')]);