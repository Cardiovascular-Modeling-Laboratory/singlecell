% 必要なデータをプロットするためのスクリプト

% ディレクトリ内の .mat ファイルを取得
dataDir = pwd; % 今の作業ディレクトリ
filePattern = fullfile(dataDir, 'fig9_v1_nuclei_stretch_later_simulation_results_*.mat');
resultFiles = dir(filePattern);

% 保存用変数
final_values = [];
mean_angles = [];
OOP_values = [];
traction_forces = [];
nucleus_obstruction = []; % 核の障害の有無を保存するための変数

% データを収���
processed_final_values = []; % 重複を避けるためのリスト
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);
    
    % ファイル名から final, nucleus_obstruction を抽出
    tokens = regexp(resultFiles(k).name, 'first(\d+\.\d+)_nuc_rel_vec_(\d+)', 'tokens');
    if ~isempty(tokens)
        final_value = str2double(tokens{1}{1});
        nuc_obstruction = str2double(tokens{1}{2});
        disp(['Final value: ', num2str(final_value), ', Nucleus Obstruction: ', num2str(nuc_obstruction)]);
    else
        disp(['Could not parse final value or nucleus obstruction from filename: ', resultFiles(k).name]);
        continue; % 次のファイルへ
    end
    
    % データ読み込み
    load(filename, 'net_res', 'F_res', 'F_store');  % 必要なデータを読み込む

    % ---- 必要なデータを計算・保存する ----
    % 例: 平均角度や OOP の計算
    mean_angle = calculate_last_mean_angle(net_res);
    OOP = calculate_OOP(net_res);
    traction_force = calculate_last_traction_force(F_store);
    
    % 結果を保存
    final_values = [final_values; final_value];
    mean_angles = [mean_angles; mean_angle];
    OOP_values = [OOP_values; OOP];
    traction_forces = [traction_forces; traction_force];
    nucleus_obstruction = [nucleus_obstruction; nuc_obstruction]; % 核の障害の有無を保存
end

% データをソート
[sorted_final_values, sort_idx] = sort(final_values);
sorted_mean_angles = mean_angles(sort_idx);
sorted_OOP_values = OOP_values(sort_idx);
sorted_traction_forces = traction_forces(sort_idx);
sorted_nucleus_obstruction = nucleus_obstruction(sort_idx);

% 核の障害の有無ごとにデータを分割
no_obstruction_idx = sorted_nucleus_obstruction == 0;
with_obstruction_idx = sorted_nucleus_obstruction == 1;

% x軸の位置を調整
x_no_obstruction = sorted_final_values(no_obstruction_idx) - 0.35;
x_with_obstruction = sorted_final_values(with_obstruction_idx) + 0.35;

% ---- グラフをプロット ----
figure;

% グラフ1: Mean Angle vs Final Value
subplot(1, 3, 1);
hold on;
bar(x_no_obstruction, sorted_mean_angles(no_obstruction_idx), 'FaceColor', 'b', 'BarWidth', 0.3);
bar(x_with_obstruction, sorted_mean_angles(with_obstruction_idx), 'FaceColor', 'r', 'BarWidth', 0.3);
xlabel('Final Value');
ylabel('Mean Angle (degrees)');
ylim([0, 90]); % 角度の範囲を 0-90 に設定
legend('No Obstruction', 'With Obstruction');
title('Mean Angle');

% グラフ2: OOP vs Final Value
subplot(1, 3, 2);
hold on;
bar(x_no_obstruction, sorted_OOP_values(no_obstruction_idx), 'FaceColor', 'b', 'BarWidth', 0.3);
bar(x_with_obstruction, sorted_OOP_values(with_obstruction_idx), 'FaceColor', 'r', 'BarWidth', 0.3);
xlabel('Final Value');
ylabel('OOP');
ylim([0, 1]); % OOPの範囲を 0-1 に設定
legend('No Obstruction', 'With Obstruction');
title('OOP');

% グラフ3: Traction Force vs Final Value
subplot(1, 3, 3);
hold on;
bar(x_no_obstruction, sorted_traction_forces(no_obstruction_idx), 'FaceColor', 'b', 'BarWidth', 0.3);
bar(x_with_obstruction, sorted_traction_forces(with_obstruction_idx), 'FaceColor', 'r', 'BarWidth', 0.3);
xlabel('Final Value');
ylabel('Traction Force (N)');
legend('No Obstruction', 'With Obstruction');
title('Traction Force');

% 全体のレイアウト調整
sgtitle('Stretch Results Analysis with and without Nucleus Obstruction');
