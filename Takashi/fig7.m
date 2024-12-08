% 必要なデータをプロットするためのスクリプト

% ディレクトリ内の .mat ファイルを取得
dataDir = pwd; % 今の作業ディレクトリ
filePattern = fullfile(dataDir, 'fig7_v1_stretch_later_simulation_results_*.mat');
resultFiles = dir(filePattern);

% 保存用変数
final_values = [];
mean_angles = [];
OOP_values = [];
traction_forces = [];

% データを収集
processed_final_values = []; % 重複を避けるためのリスト
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);
    
    % ファイル名から final を抽出
    tokens = regexp(resultFiles(k).name, 'final(\d+\.\d+)', 'tokens');
    if ~isempty(tokens)
        final_value = str2double(tokens{1}{1});
    else
        disp(['Could not parse final value from filename: ', resultFiles(k).name]);
        continue; % 次のファイルへ
    end
    
    % 重複チェック
    if ismember(final_value, processed_final_values)
        continue; % 既に処理済みの final 値はスキップ
    end
    processed_final_values = [processed_final_values, final_value];
    
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
end

% ---- グラフをプロット ----
figure;

% グラフ1: Mean Angle vs Final
subplot(1, 3, 1);
bar(final_values, mean_angles);
xlabel('Final');
ylabel('Mean Angle (degrees)');
ylim([0, 90]); % 角度の範囲を 0-90 に設定

% グラフ2: OOP vs Final
subplot(1, 3, 2);
bar(final_values, OOP_values);
xlabel('Final');
ylabel('OOP');
ylim([0, 1]); % OOPの範囲を 0-1 に設定
title('OOP');

% グラフ3: Traction Force vs Final
subplot(1, 3, 3);
bar(final_values, traction_forces);
xlabel('Final');
ylabel('Traction Force (N)');

% 全体のレイアウト調整
sgtitle('Stretch Results Analysis');
