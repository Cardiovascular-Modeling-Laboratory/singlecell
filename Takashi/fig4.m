% 必要なデータをプロットするためのスクリプト

% ディレクトリ内の .mat ファイルを取得
dataDir = pwd; % 今の作業ディレクトリ
filePattern = fullfile(dataDir, 'fig4*.mat');
resultFiles = dir(filePattern);

% 保存用変数
fold_times_str = [];
fold_times_num = [];
stretch_durations = []; % 新たに stretch duration を保存
mean_angles = [];
OOP_values = [];
traction_forces = [];

% 固定パラメータ
max_stretch_factor = 3.0; % 最大ストレッチファクター
total_time = 72 * 3600;  % 72時間分を秒に変換

% データを収集
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);
    
    % ファイル名から fold_time を抽出
    tokens = regexp(resultFiles(k).name, 'results_(\d+\.\d+)_foldtime_(\d+)', 'tokens');
    if ~isempty(tokens)
        fold_time_str = tokens{1}{2}; % fold_time を文字列として取得
        fold_time_num = str2double(fold_time_str); % fold_time を数値として取得
    else
        disp(['Could not parse fold_time from filename: ', resultFiles(k).name]);
        continue; % 次のファイルへ
    end
    
    % Stretch Duration を計算
    stretch_duration = total_time - fold_time_num; % Stretch Duration = 72*3601 - fold_time
    
    % データ読み込み
    load(filename, 'net_res', 'F_res', 'F_store');  % 必要なデータを読み込む

    % ---- 必要なデータを計算・保存する ----
    % 平均角度の計算
    mean_angle = calculate_mean_angle(net_res);
    
    % OOP の計算
    OOP = calculate_OOP(net_res);
    
    % 最大牽引力の計算
    traction_force = calculate_max_traction_force(F_store);
    
    % 結果を保存
    fold_times_str = [fold_times_str; {fold_time_str}]; % 文字列として保存
    fold_times_num = [fold_times_num; fold_time_num];   % 数値として保存
    stretch_durations = [stretch_durations; stretch_duration]; % Stretch Duration を保存
    mean_angles = [mean_angles; mean_angle];
    OOP_values = [OOP_values; OOP];
    traction_forces = [traction_forces; traction_force];
end

% ---- stretch duration でソート ----
[stretch_durations_sorted, sort_idx] = sort(stretch_durations); % Stretch Duration でソート
fold_times_str_sorted = fold_times_str(sort_idx);               % fold_time 文字列も並び替え
mean_angles_sorted = mean_angles(sort_idx);
OOP_values_sorted = OOP_values(sort_idx);
traction_forces_sorted = traction_forces(sort_idx);

% x 軸を categorical 型に変換（ソート後）
stretch_durations_categorical = categorical(string(stretch_durations_sorted)); % Stretch Duration をカテゴリに変換
stretch_durations_categorical = reordercats(stretch_durations_categorical, string(stretch_durations_sorted));

% ---- グラフをプロット ----
figure;

% グラフ1: Mean Angle vs Stretch Duration
subplot(1, 3, 1);
bar(stretch_durations_categorical, mean_angles_sorted);
xlabel('Stretch Duration (s)');
ylabel('Mean Angle (°)'); % 単位を追加
% title('Mean Angle vs Stretch Duration');
ylim([0, 90]); % Y軸を0-90°に制限

% グラフ2: OOP vs Stretch Duration
subplot(1, 3, 2);
bar(stretch_durations_categorical, OOP_values_sorted);
xlabel('Stretch Duration (s)');
ylabel('OOP (unitless)'); % 単位を追加
% title('OOP vs Stretch Duration');

% グラフ3: Traction Force vs Stretch Duration
subplot(1, 3, 3);
bar(stretch_durations_categorical, traction_forces_sorted);
xlabel('Stretch Duration (s)');
ylabel('Traction Force (N)'); % 単位を追加
% title('Traction Force vs Stretch Duration');

% 全体のレイアウト調整
sgtitle('Stretch Duration Results Analysis');


% ---- 計算用の関数群 ----

function max_force = calculate_max_traction_force(F_store)
    % 最大牽引力を計算する関数
    % F_store: 力のデータ (Nx2 の行列、[F_x, F_y])
    % max_force: 最大牽引力

    % 各時点での牽引力の大きさを計算
    forces = sqrt(F_store(:, :, 1).^2 + F_store(:, :, 2).^2); % 合力
    max_force = max(forces(:)); % 最大値を取得
end

function mean_angle = calculate_mean_angle(net_res)
    % Mean Angle を計算する関数
    % net_res: シミュレーション結果の繊維データ
    % mean_angle: 平均角度（ラジアン）

    angles = [];

    % net_res{1,1} の最終ステップのデータを使用
    net1 = net_res{1, 1};
    net = net1{end};

    for i = 1:length(net)
        net_temp = net{i};

        % net_temp の各繊維について処理
        for j = 1:size(net_temp, 3)
            % 繊維の方向ベクトルを取得
            x_start = net_temp(1, 1, j);
            y_start = net_temp(2, 1, j);
            x_end = net_temp(1, end, j);
            y_end = net_temp(2, end, j);

            % ベクトルの方向を計算
            dx = x_end - x_start;
            dy = y_end - y_start;
            angle = atan2(dy, dx); % 角度を計算
            degrees = rad2deg(angle); % ラジアンから度に変換
            angles = [angles, degrees];
        end
    end

    % 平均角度を計算
    if ~isempty(angles)
        mean_angle = mean(angles);
    else
        mean_angle = NaN;
    end
end

function OOP = calculate_OOP(net_res)
    % OOP を計算する関数
    % net_res: シミュレーション結果の繊維データ
    % OOP: Orientation Order Parameter (スカラー値)

    tensor_sum = zeros(2, 2);
    fiber_count = 0;

    % net_res{1,1} の最終ステップのデータを使用
    net1 = net_res{1, 1};
    net = net1{end};

    for i = 1:length(net)
        net_temp = net{i};
        
        % net_temp の各繊維について処理
        for j = 1:size(net_temp, 3)
            % 繊維の方向ベクトルを取得
            x_start = net_temp(1, 1, j);
            y_start = net_temp(2, 1, j);
            x_end = net_temp(1, end, j);
            y_end = net_temp(2, end, j);

            % 方向ベクトルを単位ベクトルに正規化
            dx = x_end - x_start;
            dy = y_end - y_start;
            norm_factor = sqrt(dx^2 + dy^2);
            if norm_factor == 0
                continue; % 長さゼロはスキップ
            end
            p_x = dx / norm_factor;
            p_y = dy / norm_factor;

            % テンソルを計算
            tensor = [p_x * p_x, p_x * p_y; p_x * p_y, p_y * p_y];
            tensor_sum = tensor_sum + tensor;
            fiber_count = fiber_count + 1;
        end
    end

    % 平均テンソルを計算して最大固有値を求める
    if fiber_count > 0
        mean_tensor = 2 * (tensor_sum / fiber_count) - eye(2);
        eigenvalues = eig(mean_tensor);
        OOP = max(eigenvalues);
    else
        OOP = NaN;
    end
end
