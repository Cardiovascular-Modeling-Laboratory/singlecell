% 必要なデータをプロットするためのスクリプト

% ディレクトリ内の .mat ファイルを取得
dataDir = pwd; % 今の作業ディレクトリ
filePattern = fullfile(dataDir, 'fig8_v1_stretch_later_simulation_results_*.mat');
resultFiles = dir(filePattern);

% 保存用変数
final_values = [];
mean_angles = [];
OOP_values = [];
traction_forces = [];
fold_times = []; % fold_time を保存するための変数
total_time = 72 * 3600;  % 72時間分を秒に変換

% データを収集
processed_final_values = []; % 重複を避けるためのリスト
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);
    
    % ファイル名から final と fold_time を抽出
    tokens = regexp(resultFiles(k).name, 'first(\d+\.\d+)_foldtime_(\d+)', 'tokens');
    if ~isempty(tokens)
        final_value = str2double(tokens{1}{1});
        fold_time = str2double(tokens{1}{2});
        fold_time = total_time - fold_time;
        disp(['Final value: ', num2str(final_value), ', Fold time: ', num2str(fold_time)]);
    else
        disp(['Could not parse final value or fold time from filename: ', resultFiles(k).name]);
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
    fold_times = [fold_times; fold_time]; % fold_time を保存
end

% fold_times を文字列としてソート
[sorted_fold_times, sort_idx] = sort(fold_times);
sorted_fold_times = string(sorted_fold_times);
sorted_mean_angles = mean_angles(sort_idx);
sorted_OOP_values = OOP_values(sort_idx);
sorted_traction_forces = traction_forces(sort_idx);

% ---- グラフをプロット ----
figure;

% グラフ1: Mean Angle vs Shrink Duration
subplot(1, 3, 1);
bar(sorted_fold_times, sorted_mean_angles);
xlabel('Shrink Duration');
ylabel('Mean Angle (degrees)');
ylim([0, 90]); % 角度の範囲を 0-90 に設定

% グラフ2: OOP vs Shrink Duration
subplot(1, 3, 2);
bar(sorted_fold_times, sorted_OOP_values);
xlabel('Shrink Duration');
ylabel('OOP');
ylim([0, 1]); % OOPの範囲を 0-1 に設定
title('OOP');

% グラフ3: Traction Force vs Shrink Duration
subplot(1, 3, 3);
bar(sorted_fold_times, sorted_traction_forces);
xlabel('Shrink Duration');
ylabel('Traction Force (N)');

% 全体のレイアウト調整
sgtitle('Stretch Results Analysis');
