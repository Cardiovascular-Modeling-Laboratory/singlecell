% filepath: /Users/inagakit/Documents/UCIrvine/AnnaGrosberg/singlecell/Takashi/fig3_v2.m
% 必要なデータをプロットするためのスクリプト

% ...existing code...
dataDir = pwd; % 今の作業ディレクトリ
% ディレクトリ内の .mat ファイルを取得（n=3 の結果を含むファイルパターンに変更）
filePattern = fullfile(dataDir, 'Takashi', 'results', 'fig3_v2_n=*_stretch_*.mat');
resultFiles = dir(filePattern);

resultFiles

% 保存用変数（データを蓄積するために変更）
stretch_factors = [];
mean_angles_all = [];
OOP_values_all = [];
traction_forces_all = [];

% データを収集
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);
    
    % ファイル名から stretch_factor と n を抽出
    tokens = regexp(resultFiles(k).name, 'fig3_v2_n=(\d+)_stretch_(\d+)', 'tokens');
    if ~isempty(tokens)
        n_sim = str2double(tokens{1}{1});
        stretch_factor = str2double(tokens{1}{2});
    else
        disp(['Could not parse n or stretch_factor from filename: ', resultFiles(k).name]);
        continue; % 次のファイルへ
    end
    
    % データ読み込み
    load(filename, 'net_res', 'F_res', 'F_store');  % 必要なデータを読み込む

    % ---- 必要なデータを計算・保存する ----
    % 複数のシミュレーション結果から平均角度や OOP を計算
    mean_angles = calculate_last_mean_angle(net_res);      % サイズ: [n_sim x 1]
    OOP_values = calculate_OOP(net_res);        % サイズ: [n_sim x 1]
    traction_forces = calculate_last_traction_force(F_store); % サイズ: [n_sim x 1]
    
    % 結果を保存
    stretch_factors = [stretch_factors; stretch_factor];
    mean_angles_all = [mean_angles_all; mean_angles];
    OOP_values_all = [OOP_values_all; OOP_values];
    traction_forces_all = [traction_forces_all; traction_forces];
end

% ユニークなストレッチファクターを取得
unique_stretch_factors = unique(stretch_factors);

% ---- データ収集部分のn数制限を追加 ----
filtered_stretch_factors = [];
filtered_mean_angles_all = [];
filtered_OOP_values_all = [];
filtered_traction_forces_all = [];

for i = 1:length(unique_stretch_factors)
    sf = unique_stretch_factors(i);
    idx = find(stretch_factors == sf); % 条件を満たすインデックスを取得
    
    % n数を3つに制限（最初の3つを選択、必要に応じてランダムサンプリングも可能）
    if length(idx) > 3
        idx = idx(1:3); % 最初の3つを選択
        % idx = idx(randperm(length(idx), 3)); % ランダムに3つ選択する場合はこちらを使用
    end
    
    % 制限後のデータを保存
    filtered_stretch_factors = [filtered_stretch_factors; stretch_factors(idx)];
    filtered_mean_angles_all = [filtered_mean_angles_all; mean_angles_all(idx)];
    filtered_OOP_values_all = [filtered_OOP_values_all; OOP_values_all(idx)];
    filtered_traction_forces_all = [filtered_traction_forces_all; traction_forces_all(idx)];
    
    % n数の表示（制限後は必ず3になる）
    n_count = length(idx); % 制限後のn数
    disp(['Stretch factor: ', num2str(sf), ', n: ', num2str(n_count)]);
end

% ---- 平均と標準偏差を計算（制限後） ----
mean_angles_mean = zeros(length(unique_stretch_factors),1);
mean_angles_sd = zeros(length(unique_stretch_factors),1);
OOP_values_mean = zeros(length(unique_stretch_factors),1);
OOP_values_sd = zeros(length(unique_stretch_factors),1);
traction_forces_mean = zeros(length(unique_stretch_factors),1);
traction_forces_sd = zeros(length(unique_stretch_factors),1);

for i = 1:length(unique_stretch_factors)
    sf = unique_stretch_factors(i);
    idx = filtered_stretch_factors == sf; % フィルタリング後のデータに対して条件を適用
    
    mean_angles_mean(i) = mean(filtered_mean_angles_all(idx));
    mean_angles_sd(i) = std(filtered_mean_angles_all(idx));
    OOP_values_mean(i) = mean(filtered_OOP_values_all(idx));
    OOP_values_sd(i) = std(filtered_OOP_values_all(idx));
    traction_forces_mean(i) = mean(filtered_traction_forces_all(idx));
    traction_forces_sd(i) = std(filtered_traction_forces_all(idx));
end



% ---- データ収集部分のn数を表示 ----
for i = 1:length(unique_stretch_factors)
    sf = unique_stretch_factors(i);
    idx = stretch_factors == sf;
    n_count = sum(idx); % 該当条件のデータ数
    disp(['Stretch factor: ', num2str(sf), ', n: ', num2str(n_count)]);
end


% Mean Angleに対するANOVAとポストホックテスト
[p_angle, tbl_angle, stats_angle] = anova1(filtered_mean_angles_all, filtered_stretch_factors, 'off');
[c_angle, m_angle, h_angle, gnames_angle] = multcompare(stats_angle, 'Display', 'off');

% OOPに対するANOVAとポストホックテスト
[p_OOP, tbl_OOP, stats_OOP] = anova1(filtered_OOP_values_all, filtered_stretch_factors, 'off');
[c_OOP, m_OOP, h_OOP, gnames_OOP] = multcompare(stats_OOP, 'Display', 'off');

% Traction Forceに対するANOVAとポストホックテスト
[p_traction, tbl_traction, stats_traction] = anova1(filtered_traction_forces_all, filtered_stretch_factors, 'off');
[c_traction, m_traction, h_traction, gnames_traction] = multcompare(stats_traction, 'Display', 'off');


debug_show_asterisks = false;
% ---- グラフをプロット ----
figure;

% グラフ1: Mean Angle vs Stretch Factor
subplot(1, 3, 1);
bar(unique_stretch_factors, mean_angles_mean);
hold on;lo
errorbar(unique_stretch_factors, mean_angles_mean, mean_angles_sd, 'LineStyle', 'none', 'LineWidth', 2); % 太さを変更
xlabel('Stretch Factor');
ylabel('Mean Angle (degrees)');
ylim([0, 360]); % 角度の範囲を 0-180 に設定

sig_positions = max(mean_angles_mean + mean_angles_sd) + 5;
significant_pairs = c_angle(c_angle(:,6) < 0.05, 1:2); % p<0.05のペア
offset = 0;
for i = 1:size(significant_pairs, 1)
    group1 = significant_pairs(i,1);
    group2 = significant_pairs(i,2);
    x1 = unique_stretch_factors(group1);
    x2 = unique_stretch_factors(group2);
    y = sig_positions + offset;
    plot([x1, x1, x2, x2], [y, y+1, y+1, y], '-k', 'LineWidth', 1);
    text(mean([x1, x2]), y+1.5, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    offset = offset + 5; % 線が重ならないように高さを調整
end

% グラフ2: OOP vs Stretch Factor
subplot(1, 3, 2);
bar(unique_stretch_factors, OOP_values_mean);
hold on;
errorbar(unique_stretch_factors, OOP_values_mean, OOP_values_sd, 'LineStyle', 'none', 'LineWidth', 2); % 太さを変更
xlabel('Stretch Factor');
ylabel('OOP');
title('OOP');
ylim([0, 1]); % OOPの範囲を 0-1 に設定

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
    offset = offset + 0.05; % 線が重ならないように高さを調整
end

% グラフ3: Traction Force vs Stretch Factor
subplot(1, 3, 3);
bar(unique_stretch_factors, traction_forces_mean);
hold on;
errorbar(unique_stretch_factors, traction_forces_mean, traction_forces_sd, 'LineStyle', 'none', 'LineWidth', 2); % 太さを変更
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
    offset = offset + 0.05; % 線が重ならないように高さを調整
end

% 全体のレイアウト調整
sgtitle('Stretch Results Analysis');

% ---- ANOVAによる統計検定 ----
% 各変数に対して一元配置分散分析を実施
[p_angle, tbl_angle] = anova1(mean_angles_all, stretch_factors, 'off');
[p_OOP, tbl_OOP] = anova1(OOP_values_all, stretch_factors, 'off');
[p_traction, tbl_traction] = anova1(traction_forces_all, stretch_factors, 'off');

% 結果を表示
disp('ANOVA for Mean Angles:');
disp(tbl_angle);
disp(['p-value: ', num2str(p_angle)]);

disp('ANOVA for OOP:');
disp(tbl_OOP);
disp(['p-value: ', num2str(p_OOP)]);

disp('ANOVA for Traction Forces:');
disp(tbl_traction);
disp(['p-value: ', num2str(p_traction)]);



