% plot_OOP_2D_heatmap_with_tensor_method.m
% このスクリプトは、指定されたディレクトリ内の各 stretch_factor と fold_time に対して
% 最新のファイルを使用して Orientation Order Parameter (OOP) をテンソルの最大固有値を用いて計算し、結果を heatmap として表示します。

% データディレクトリ
dataDir = pwd; % .matファイルが保存されているディレクトリを指定
filePattern = fullfile(dataDir, 'stretch_first_simulation_results_*.mat');
resultFiles = dir(filePattern);

% 各 stretch_factor と fold_time ごとの最新ファイルを選択
latestFiles = containers.Map();
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    
    % ファイル名から stretch_factor と fold_time を抽出
    tokens = regexp(resultFiles(k).name, 'results_(\d+\.\d+)_foldtime_(\d+)', 'tokens');
    if ~isempty(tokens)
        stretch_factor = str2double(tokens{1}{1});
        fold_time = str2double(tokens{1}{2});
        key = sprintf('%0.1f_%d', stretch_factor, fold_time); % stretch_factor と fold_time を組み合わせたキー
        
        % 最新ファイルの更新チェック
        if isKey(latestFiles, key)
            existingFile = latestFiles(key);
            if resultFiles(k).datenum > existingFile.datenum
                latestFiles(key) = resultFiles(k);
            end
        else
            latestFiles(key) = resultFiles(k);
        end
    else
        disp(['Could not parse stretch_factor or fold_time from filename: ', resultFiles(k).name]);
    end
end

% OOP 計算のためのリスト
stretch_factors = [];
fold_times = [];
OOP_values = [];

% 各最新ファイルをループして OOP を計算
latestKeys = keys(latestFiles);
for k = 1:length(latestKeys)
    key = latestKeys{k};
    resultFile = latestFiles(key);
    filename = fullfile(resultFile.folder, resultFile.name);
    disp(['Processing file: ', filename]);
    
    % キーから stretch_factor と fold_time を分解
    tokens = regexp(key, '(\d+\.\d+)_(\d+)', 'tokens');
    stretch_factor = str2double(tokens{1}{1});
    fold_time = str2double(tokens{1}{2});
    
    % データ読み込み
    load(filename, 'net_res');  % net_res が fiber データを格納していると仮定
    
    % テンソル計算のための初期化
    tensor_sum = zeros(2, 2); % テンソルの和
    fiber_count = 0;
    
    % 各 fiber のテンソル要素を計算
    net1 = net_res{1, 1};  % 各ファイルの最初のシミュレーション結果を使用
    net = net1{end}; % 最終ステップの結果
    for i = 1:length(net)
        net_temp = net{i};
        
        % net_temp の形状を確認して処理
        if size(net_temp, 1) < 2 || size(net_temp, 2) < 2
            disp('Skipping fiber due to unexpected shape.');
            continue; % 十分なデータがない場合はスキップ
        end
        
        for j = 1:size(net_temp, 3)
            % fiber の方向ベクトルを取得
            x_start = net_temp(1, 1, j);
            y_start = net_temp(2, 1, j);
            x_end = net_temp(1, end, j);
            y_end = net_temp(2, end, j);

            % 方向ベクトルを単位ベクトルに正規化
            dx = x_end - x_start;
            dy = y_end - y_start;
            norm_factor = sqrt(dx^2 + dy^2);
            if norm_factor == 0
                continue; % 長さがゼロのベクトルはスキップ
            end
            p_x = dx / norm_factor;
            p_y = dy / norm_factor;

            % テンソルを計算して和に加える
            tensor = [p_x * p_x, p_x * p_y; p_x * p_y, p_y * p_y];
            tensor_sum = tensor_sum + tensor;
            fiber_count = fiber_count + 1;
        end
    end
    
    % 平均テンソルを計算し、単位行列の影響を引く
    if fiber_count > 0
        mean_tensor = 2 * (tensor_sum / fiber_count) - eye(2);
    else
        mean_tensor = zeros(2, 2);
    end

    % 最大固有値を OOP として計算
    eigenvalues = eig(mean_tensor);
    OOP = max(eigenvalues);
    
    % OOP と対応する stretch_factor と fold_time を保存
    stretch_factors = [stretch_factors; stretch_factor];
    fold_times = [fold_times; fold_time];
    OOP_values = [OOP_values; OOP];
end

% stretch_factor と fold_time のユニークな値
unique_stretch_factors = unique(stretch_factors);
unique_fold_times = unique(fold_times);

% ヒートマップ用の OOP マトリックスを作成
OOP_matrix = NaN(length(unique_fold_times), length(unique_stretch_factors));
for i = 1:length(stretch_factors)
    sf_idx = find(unique_stretch_factors == stretch_factors(i));
    ft_idx = find(unique_fold_times == fold_times(i));
    OOP_matrix(ft_idx, sf_idx) = OOP_values(i);
end

% ヒートマップのプロット
figure;
heatmap(unique_stretch_factors, unique_fold_times, OOP_matrix, 'Colormap', jet, 'ColorLimits', [0, 1]);
xlabel('Stretch Factor');
ylabel('Fold Time (s)');
title('Orientation Order Parameter (OOP) by Stretch Factor and Fold Time (Tensor Method)');
