% calculate_OOP_per_file_fixed.m
% このスクリプトは、指定されたディレクトリ内のすべての結果ファイルに対して
% OOP を計算し、stretch_factor と fold_time に対する OOP をヒートマップとして表示します。

% データディレクトリ
dataDir = pwd; % .matファイルが保存されているディレクトリを指定
filePattern = fullfile(dataDir, 'stretch_first_simulation_results_*.mat');
resultFiles = dir(filePattern);

% OOP 計算のためのリスト
stretch_factors = [];
fold_times = [];
OOP_values = [];

% すべての結果ファイルに対してループ
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    disp(['Processing file: ', filename]);
    
    % ファイル名から stretch_factor と fold_time を抽出
    tokens = regexp(resultFiles(k).name, 'results_(\d+\.\d+)_foldtime_(\d+)', 'tokens');
    if ~isempty(tokens)
        stretch_factor = str2double(tokens{1}{1});
        fold_time = str2double(tokens{1}{2});
    else
        disp(['Could not parse stretch_factor or fold_time from filename: ', resultFiles(k).name]);
        continue; % 次のファイルへ
    end
    
    % データ読み込み
    load(filename, 'net_res');  % net_res が fiber データを格納していると仮定
    
    % シミュレーション条件の設定
    j2 = 1; % nucleus placement index
    nsim = size(net_res, 2);
    
    % 繊維セグメントの情報を格納するリストを初期化
    p_vectors = [];
    
    % 各シミュレーションに対してループ
    for m = 1:nsim
        % 繊維ネットワークのデータを取得
        net1 = net_res{j2, m};
        net = net1{end}; % 最終結果のみを使用
    
        % グリッドの設定（10x10 のグリッド）
        grid_size = 10;
        % データの範囲を取得（繊維の全座標から最小値と最大値を計算）
        all_x = [];
        all_y = [];
        for j = 1:length(net)
            net_temp = net{j};
            if isempty(net_temp)
                continue;
            end
            x_coords = net_temp(1, :, :);
            x_coords = x_coords(:); % ベクトルに変換
            all_x = [all_x; x_coords];
            y_coords = net_temp(2, :, :);
            y_coords = y_coords(:); % ベクトルに変換
            all_y = [all_y; y_coords];
        end
        if isempty(all_x) || isempty(all_y)
            disp('No fiber data found in this simulation.');
            continue; % 次のシミュレーションへ
        end
        x_min = min(all_x);
        x_max = max(all_x);
        y_min = min(all_y);
        y_max = max(all_y);
    
        % グリッドエッジを作成
        x_edges = linspace(x_min, x_max, grid_size+1);
        y_edges = linspace(y_min, y_max, grid_size+1);
    
        % 各グリッドセルについてループ
        for xi = 1:grid_size
            for yi = 1:grid_size
                % グリッドセルの範囲を定義
                x_start = x_edges(xi);
                x_end = x_edges(xi+1);
                y_start = y_edges(yi);
                y_end = y_edges(yi+1);
    
                % グリッドセル内のセグメントの方向ベクトルを格納するリスト
                p_vectors_cell = [];
    
                % 繊維ネットワークをループ
                for j = 1:length(net)
                    net_temp = net{j};
                    if isempty(net_temp)
                        continue;
                    end
                    for idx = 1:size(net_temp, 3)
                        % 繊維の座標を取得
                        net_x_temp = net_temp(1, :, idx);
                        net_y_temp = net_temp(2, :, idx);
    
                        % 各セグメント（隣接するポイント間）についてループ
                        for s = 1:length(net_x_temp)-1
                            x_seg_start = net_x_temp(s);
                            y_seg_start = net_y_temp(s);
                            x_seg_end = net_x_temp(s+1);
                            y_seg_end = net_y_temp(s+1);
    
                            % セグメントの中点を計算
                            x_mid = (x_seg_start + x_seg_end) / 2;
                            y_mid = (y_seg_start + y_seg_end) / 2;
    
                            % 中点がグリッドセル内にあるかを判定
                            if x_mid >= x_start && x_mid < x_end && y_mid >= y_start && y_mid < y_end
                                % セグメントの方向ベクトルを計算
                                dx = x_seg_end - x_seg_start;
                                dy = y_seg_end - y_seg_start;
    
                                % ベクトルの長さがゼロの場合はスキップ
                                if dx == 0 && dy == 0
                                    continue;
                                end
    
                                % ベクトルを正規化
                                magnitude = sqrt(dx^2 + dy^2);
                                p_i = [dx; dy] / magnitude;
    
                                % グリッドセル内のベクトルリストに追加
                                p_vectors_cell = [p_vectors_cell, p_i];
                            end
                        end
                    end
                end
    
                % グリッドセル内のベクトルを全体のリストに追加
                p_vectors = [p_vectors, p_vectors_cell];
            end
        end
    end
    
    % 全グリッドセルのベクトルを用いて OOP を計算
    N = size(p_vectors, 2);
    if N > 0
        % テンソルの合計を計算
        T_sum = zeros(2, 2);
        for n = 1:N
            p_i = p_vectors(:, n);
            T_i = 2 * (p_i * p_i') - eye(2);
            T_sum = T_sum + T_i;
        end
    
        % 平均テンソルを計算
        T_mean = T_sum / N;
    
        % 最大固有値を OOP として取得
        eigenvalues = eig(T_mean);
        OOP = max(eigenvalues);
    else
        OOP = NaN; % データがない場合は NaN
    end
    
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
for i = 1:length(OOP_values)
    sf_idx = find(unique_stretch_factors == stretch_factors(i));
    ft_idx = find(unique_fold_times == fold_times(i));
    OOP_matrix(ft_idx, sf_idx) = OOP_values(i);
end

% ヒートマップのプロット
figure;
h = heatmap(unique_stretch_factors, unique_fold_times, OOP_matrix, 'Colormap', jet, 'ColorLimits', [0, 1]);
xlabel('Stretch Factor');
ylabel('Fold Time (s)');
title('Orientation Order Parameter (OOP) by Stretch Factor and Fold Time');
h.ColorbarVisible = 'on';
