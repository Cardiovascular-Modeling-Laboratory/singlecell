function OOP = calculate_OOP(net_res)
    % OOP を計算する関数
    % net_res: シミュレーション結果の繊維データ
    % OOP: Orientation Order Parameter (スカラー値)

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
    OOP = calculate_OOP_from_vectors(total_vectors);
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
function OOP = calculate_OOP_from_vectors(vectors)
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
