% 必要なデータをプロットするためのスクリプト

% ディレクトリ内の .mat ファイルを取得
dataDir = pwd; % 今の作業ディレクトリ
filePattern = fullfile(dataDir, 'fig3_v2*.mat');
resultFiles = dir(filePattern);

% 保存用変数
stretch_factors = [];
mean_angles = [];
OOP_values = [];
traction_forces = [];

% データを収集
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);
    
    % ファイル名から stretch_factor を抽出
    tokens = regexp(resultFiles(k).name, 'results_(\d+\.\d+)_foldtime_(\d+)', 'tokens');
    if ~isempty(tokens)
        stretch_factor = str2double(tokens{1}{1});
        fold_time = str2double(tokens{1}{2});
    else
        disp(['Could not parse stretch_factor or fold_time from filename: ', resultFiles(k).name]);
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
    stretch_factors = [stretch_factors; stretch_factor];
    mean_angles = [mean_angles; mean_angle];
    OOP_values = [OOP_values; OOP];
    traction_forces = [traction_forces; traction_force];
end

% ---- グラフをプロット ----
figure;

% グラフ1: Mean Angle vs Stretch Length
subplot(1, 3, 1);
bar(stretch_factors, mean_angles);
xlabel('Stretch Factor');
ylabel('Mean Angle (degrees)');
ylim([0, 90]); % 角度の範囲を 0-180 に設定

% グラフ2: OOP vs Stretch Factor
subplot(1, 3, 2);
bar(stretch_factors, OOP_values);
xlabel('Stretch Factor');
ylabel('OOP');
title('OOP');
ylim([0, 1]); % OOPの範囲を 0-1 に設定

% グラフ3: Traction Force vs Stretch Factor
subplot(1, 3, 3);
bar(stretch_factors, traction_forces);
xlabel('Stretch Factor');
ylabel('Traction Force (N)');

% 全体のレイアウト調整
sgtitle('Stretch Results Analysis');

