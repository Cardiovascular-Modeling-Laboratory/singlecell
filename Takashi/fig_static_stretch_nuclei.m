% filepath: /Users/inagakit/Documents/UCIrvine/AnnaGrosberg/singlecell/Takashi/fig9_v2.m
% 必要なデータをプロットするためのスクリプト

% ディレクトリ内の .mat ファイルを取得（n=3 の結果を含むファイルパターンに変更）
dataDir = pwd;
filePattern = fullfile(dataDir, 'Takashi', 'results', 'fig9_v2_n=*_final*.mat');
resultFiles = dir(filePattern);

% 保存用変数（データを蓄積するために変更）
final_values = [];
mean_angles_all = [];
OOP_values_all = [];
traction_forces_all = [];
nucleus_obstruction_all = [];

% データを収集
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);

    % ファイル名から n、final_value、nucleus_obstruction を抽出
    tokens = regexp(resultFiles(k).name, 'fig9_v2_n=(\d+)_final(\d+)_nuc_rel_vec_(\d+)', 'tokens');
    if ~isempty(tokens)
        n_sim = str2double(tokens{1}{1});
        final_value = str2double(tokens{1}{2});
        nuc_obstruction = str2double(tokens{1}{3});
    else
        disp(['Could not parse n, final value, or nucleus obstruction from filename: ', resultFiles(k).name]);
        continue; % 次のファイルへ
    end

    % データ読み込み
    load(filename, 'net_res', 'F_res', 'F_store');  % 必要なデータを読み込む

    % ---- 必要なデータを計算・保存する ----
    % 複数のシミュレーション結果から平均角度や OOP を計算
    mean_angles = calculate_last_mean_angle(net_res);      % サイズ: [n_sim x 1]
    OOP_values = calculate_OOP(net_res);                   % サイズ: [n_sim x 1]
    traction_forces = calculate_last_traction_force(F_store); % サイズ: [n_sim x 1]

    % 結果を保存
    final_values = [final_values; final_value];
    mean_angles_all = [mean_angles_all; mean_angles];
    OOP_values_all = [OOP_values_all; OOP_values];
    traction_forces_all = [traction_forces_all; traction_forces];
    nucleus_obstruction_all = [nucleus_obstruction_all; nuc_obstruction];
end

% データをソート
[sorted_final_values, sort_idx] = sort(final_values);
sorted_mean_angles = mean_angles_all(sort_idx);
sorted_OOP_values = OOP_values_all(sort_idx);
sorted_traction_forces = traction_forces_all(sort_idx);
sorted_nucleus_obstruction = nucleus_obstruction_all(sort_idx);

% 核の障害の有無ごとにデータを分割
no_obstruction_idx = sorted_nucleus_obstruction == 0;
with_obstruction_idx = sorted_nucleus_obstruction == 1;

% ユニークなファイナル値を取得
unique_final_values = unique(sorted_final_values);

% 各ファイナル値ごとにデータを分け、核の有無で平均と標準偏差を計算
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

% x軸の位置を調整
x_no_obstruction = unique_final_values;
x_with_obstruction = unique_final_values + 0.35;

% デバッグ用のフラグ
debug_show_asterisks = true;

% ---- グラフをプロット ----
figure;

% グラフ1: Mean Angle vs Final Value
subplot(1, 3, 1);
hold on;
bar1 = bar(x_no_obstruction, mean_angles_mean_no, 'FaceColor', 'b', 'BarWidth', 0.3);
errorbar(x_no_obstruction, mean_angles_mean_no, mean_angles_sd_no, '.', 'Color', 'k');
bar2 = bar(x_with_obstruction, mean_angles_mean_with, 'FaceColor', 'r', 'BarWidth', 0.3);
errorbar(x_with_obstruction, mean_angles_mean_with, mean_angles_sd_with, '.', 'Color', 'k');
xlabel('Final Value');
ylabel('Mean Angle (degrees)');
ylim([0, 90]); % 角度の範囲を 0-90 に設定
legend([bar1, bar2], {'No Obstruction', 'With Obstruction'});
title('Mean Angle');

% ANOVAの結果に基づいて有意差を表示
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

% グラフ2: OOP vs Final Value
subplot(1, 3, 2);
hold on;
bar1 = bar(x_no_obstruction, OOP_values_mean_no, 'FaceColor', 'b', 'BarWidth', 0.3);
errorbar(x_no_obstruction, OOP_values_mean_no, OOP_values_sd_no, '.', 'Color', 'k');
bar2 = bar(x_with_obstruction, OOP_values_mean_with, 'FaceColor', 'r', 'BarWidth', 0.3);
errorbar(x_with_obstruction, OOP_values_mean_with, OOP_values_sd_with, '.', 'Color', 'k');
xlabel('Final Value');
ylabel('OOP');
ylim([0, 1]); % OOPの範囲を 0-1 に設定
legend([bar1, bar2], {'No Obstruction', 'With Obstruction'});
title('OOP');

% ANOVAの結果に基づいて有意差を表示
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

% グラフ3: Traction Force vs Final Value
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

% ANOVAの結果に基づいて有意差を表示
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

% 全体のレイアウト調整
sgtitle('Stretch Results Analysis with and without Nucleus Obstruction');

% ---- ANOVAによる統計検定 ----
% 二元配置分散分析を実施
group = {sorted_final_values, sorted_nucleus_obstruction};

% Mean Angleに対するANOVA
[p_angle, tbl_angle, stats_angle] = anovan(sorted_mean_angles, group, 'model', 'interaction', 'varnames', {'FinalValue', 'NucleusObstruction'});

% OOPに対するANOVA
[p_OOP, tbl_OOP, stats_OOP] = anovan(sorted_OOP_values, group, 'model', 'interaction', 'varnames', {'FinalValue', 'NucleusObstruction'});

% Traction Forceに対するANOVA
[p_traction, tbl_traction, stats_traction] = anovan(sorted_traction_forces, group, 'model', 'interaction', 'varnames', {'FinalValue', 'NucleusObstruction'});

% 結果を表示
disp('ANOVA for Mean Angles:');
disp(tbl_angle);
disp(['p-values: ', num2str(p_angle')]);

disp('ANOVA for OOP:');
disp(tbl_OOP);
disp(['p-values: ', num2str(p_OOP')]);

disp('ANOVA for Traction Forces:');
disp(tbl_traction);
disp(['p-values: ', num2str(p_traction')]);