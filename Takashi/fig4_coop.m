% メインスクリプト

% データディレクトリ
dataDir = pwd; % 現在のディレクトリを使用
filePattern = fullfile(dataDir, 'fig4_v1_*.mat');
resultFiles = dir(filePattern);
total_time = 72 * 3600;  % 72時間分を秒に変換
% ファイル情報の抽出
file_info = struct('filename', {}, 'fold_time', {}, 'stretch_length', {});

for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    % ファイル名からfold_timeとstretch_lengthを抽出
    tokens = regexp(resultFiles(k).name, 'stretch_later_simulation_results_(\d+\.\d+)_foldtime_(\d+)', 'tokens');
    if ~isempty(tokens)
        stretch_length = str2double(tokens{1}{1});
        fold_time = str2double(tokens{1}{2});
        fold_time = total_time - fold_time;
        file_info(end+1) = struct('filename', filename, 'fold_time', fold_time, 'stretch_length', stretch_length);
    else
        disp(['Could not parse stretch_length or fold_time from filename: ', filename]);
    end
end

% ファイルペアの生成
n_files = length(file_info);
COOP_matrix = NaN(n_files, n_files); % fold_timeごとのヒートマップ

for i = 1:n_files
    for j = i:n_files
        if i == j, continue; end % 同じfold_timeの組み合わせはスキップ
        
        % ファイルペアの読み込み
        file1 = file_info(i).filename;
        file2 = file_info(j).filename;
        [OOP1, grid1] = calculate_OOP_and_fiber_orientation_grid(file1);
        [OOP2, grid2] = calculate_OOP_and_fiber_orientation_grid(file2);
        
        % COOPの計算
        COOP_matrix(i, j) = calculate_COOP(OOP1, OOP2, grid1, grid2);
    end
end

% fold_timesを昇順にソート
fold_times = [file_info.fold_time];
[sorted_fold_times, sort_idx] = sort(fold_times);

% COOPマトリックスを並び替え
sorted_COOP_matrix = COOP_matrix(sort_idx, sort_idx);

% 対角成分以外のすべての成分が埋まっているか確認しつつ、対称化
n = size(sorted_COOP_matrix, 1);
for i = 1:n
    for j = i+1:n
        if isnan(sorted_COOP_matrix(i, j)) && isnan(sorted_COOP_matrix(j, i))
            error('Error: NaN found in off-diagonal elements at (%d, %d)', i, j);
        elseif isnan(sorted_COOP_matrix(i, j)) % i, j の要素が NaN の場合
            sorted_COOP_matrix(i, j) = sorted_COOP_matrix(j, i);
        elseif isnan(sorted_COOP_matrix(j, i)) % j, i の要素が NaN の場合
            sorted_COOP_matrix(j, i) = sorted_COOP_matrix(i, j);
        end
    end
end

% 再確認：対角成分以外にNaNが残っていないかチェック
off_diag_elements = sorted_COOP_matrix(~eye(size(sorted_COOP_matrix)));
if any(isnan(off_diag_elements))
    error('Error: There are NaN values in the off-diagonal elements.');
end

% 完全な行列を表示
disp('Completed COOP matrix:');
disp(sorted_COOP_matrix);

% ヒートマップのプロット（ソート済みfold_timesを使用）
figure;
h = heatmap(sorted_fold_times, sorted_fold_times, sorted_COOP_matrix, ...
    'Colormap', jet, 'ColorLimits', [0, 1]);

xlabel('Fold Time 1 (s)');
ylabel('Fold Time 2 (s)');
title('COOP Heatmap (Sorted Fold Times)');
h.ColorbarVisible = 'on';




% 結果をテキストで保存
T = array2table(COOP_matrix, 'RowNames', string(fold_times), 'VariableNames', string(fold_times));
writetable(T, 'COOP_results.csv', 'WriteRowNames', true);

%% 関数: OOPとグリッドごとの繊維方向を計算
function [OOP, grid] = calculate_OOP_and_fiber_orientation_grid(filename)
    load(filename, 'net_res');
    grid_size = 10; % グリッドの数
    grid = zeros(grid_size, grid_size, 2); % 平均方向ベクトルを格納
    total_vectors = []; % 全繊維方向ベクトル
    
    for sim_idx = 1:length(net_res)
        net = net_res{1, sim_idx}{end}; % 最終ステップ
        % グリッドごとにベクトルを収集
        [grid_vectors, vectors] = collect_grid_vectors(net, grid_size);
        grid = grid + grid_vectors;
        total_vectors = [total_vectors, vectors];
    end
    
    % OOPの計算
    OOP = calculate_OOP(total_vectors);
end

%% 関数: COOPを計算
function COOP = calculate_COOP(OOP1, OOP2, grid1, grid2)
    N = size(grid1, 1) * size(grid1, 2); % グリッド数
    T_sum = zeros(2, 2);
    
    for xi = 1:size(grid1, 1)
        for yi = 1:size(grid1, 2)
            % グリッドベクトルを取得
            v1 = squeeze(grid1(xi, yi, :));
            v2 = squeeze(grid2(xi, yi, :));
            
            if any(v1) && any(v2) % 非ゼロベクトルのみ計算
                % ベクトルを3次元に拡張（z成分を0とする）
                v1_3d = [v1; 0];
                v2_3d = [v2; 0];
                
                % cos(θ)とsin(θ)を計算
                f_x = dot(v1, v2); % cos(θ)
                f_y = norm(cross(v1_3d, v2_3d)); % sin(θ)
                
                % f_iを計算
                f_i = [f_x; f_y];
                
                % T_iを計算
                T_i = 2 * (f_i * f_i') - eye(2);
                T_sum = T_sum + T_i;
            end
        end
    end
    
    % 平均テンソルを計算
    T_mean = T_sum / N;
    COOP = max(eig(T_mean)); % 最大固有値
end

%% 関数: グリッドごとの繊維方向ベクトルを収集
function [grid_vectors, total_vectors] = collect_grid_vectors(net, grid_size)
    % グリッドエッジの初期化
    all_x = [];
    all_y = [];
    for j = 1:length(net)
        net_temp = net{j};
        if isempty(net_temp)
            continue;
        end
        % 繊維データが3次元配列の場合を考慮
        if size(net_temp, 1) >= 2
            x_coords = net_temp(1, :, :);
            y_coords = net_temp(2, :, :);
            % 2次元配列に変換（多次元配列をベクトルに整形）
            x_coords = x_coords(:);
            y_coords = y_coords(:);
            all_x = [all_x; x_coords];
            all_y = [all_y; y_coords];
        else
            warning('Unexpected net_temp dimensions in net{%d}', j);
        end
    end

    % グリッドエッジを定義
    x_min = min(all_x); x_max = max(all_x);
    y_min = min(all_y); y_max = max(all_y);
    x_edges = linspace(x_min, x_max, grid_size+1);
    y_edges = linspace(y_min, y_max, grid_size+1);

    % グリッドごとのベクトル収集
    grid_vectors = zeros(grid_size, grid_size, 2); % 平均ベクトル格納
    total_vectors = [];
    
    for xi = 1:grid_size
        for yi = 1:grid_size
            % グリッドセルの範囲
            x_start = x_edges(xi); x_end = x_edges(xi+1);
            y_start = y_edges(yi); y_end = y_edges(yi+1);
            
            % グリッドセル内のベクトル
            cell_vectors = [];
            
            for j = 1:length(net)
                net_temp = net{j};
                if isempty(net_temp)
                    continue;
                end
                
                for idx = 1:size(net_temp, 2)-1
                    % 繊維の始点と終点
                    if size(net_temp, 1) >= 2
                        p1 = net_temp(:, idx);
                        p2 = net_temp(:, idx+1);
                        mid_point = (p1 + p2) / 2;

                        % 中点がグリッド内にあるか判定
                        if mid_point(1) >= x_start && mid_point(1) < x_end && ...
                           mid_point(2) >= y_start && mid_point(2) < y_end
                            % ベクトルを計算
                            vec = p2 - p1;
                            vec = vec / norm(vec); % 正規化
                            cell_vectors = [cell_vectors, vec];
                        end
                    else
                        warning('Unexpected net_temp dimensions in net{%d}', j);
                    end
                end
            end
            
            % グリッド平均ベクトルを計算
            if ~isempty(cell_vectors)
                avg_vec = mean(cell_vectors, 2);
                avg_vec = avg_vec / norm(avg_vec);
                grid_vectors(xi, yi, :) = avg_vec;
                total_vectors = [total_vectors, avg_vec];
            end
        end
    end
end


%% 関数: OOPを計算
function OOP = calculate_OOP(vectors)
    N = size(vectors, 2);
    if N > 0
        T_sum = zeros(2, 2);
        for i = 1:N
            v = vectors(:, i);
            T_i = 2 * (v * v') - eye(2);
            T_sum = T_sum + T_i;
        end
        T_mean = T_sum / N;
        OOP = max(eig(T_mean));
    else
        OOP = NaN;
    end
end
