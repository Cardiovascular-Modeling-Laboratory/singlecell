% filepath: /Users/inagakit/Documents/UCIrvine/AnnaGrosberg/singlecell/Takashi/fig9_coop.m

% fig9_coop.m
% このスクリプトは、核あり・なしのデータセットペア間でCOOPを計算し、バーグラフで結果をプロットします。

%% データディレクトリとファイルパターンの設定
dataDir = pwd;
filePattern = fullfile(dataDir, 'Takashi', 'results', 'fig9_v2_n=*_final*_nuc_rel_vec_*.mat');
resultFiles = dir(filePattern);

% ファイル情報の初期化
file_info = struct('filename', {}, 'n_sim', {}, 'final_value', {}, 'nucleus_obstruction', {});

% ファイル情報の収集
for k = 1:length(resultFiles)
    filename = fullfile(resultFiles(k).folder, resultFiles(k).name);
    % ファイル名から必要な情報を抽出
    tokens = regexp(resultFiles(k).name, 'fig9_v2_n=(\d+)_final(\d+)_nuc_rel_vec_(\d+)', 'tokens');
    if ~isempty(tokens)
        n_sim = str2double(tokens{1}{1});
        final_value = str2double(tokens{1}{2});
        nucleus_obstruction = str2double(tokens{1}{3});
        file_info(end+1) = struct('filename', filename, 'n_sim', n_sim, 'final_value', final_value, 'nucleus_obstruction', nucleus_obstruction);
    else
        disp(['Could not parse filename: ', resultFiles(k).name]);
        continue;
    end
end

% ユニークなストレッチ条件を取得
final_values = [file_info.final_value];
unique_final_values = unique(final_values);

% COOP結果の初期化
COOP_results = [];

% 各ストレッチ条件ごとに処理
for i = 1:length(unique_final_values)
    fv = unique_final_values(i);
    % 核ありとなしのファイルを取得
    idx_no_nucleus = find([file_info.final_value] == fv & [file_info.nucleus_obstruction] == 0);
    idx_with_nucleus = find([file_info.final_value] == fv & [file_info.nucleus_obstruction] == 1);

    % ペアの数を決定
    num_pairs = min(length(idx_no_nucleus), length(idx_with_nucleus));
    for j = 1:num_pairs
        % 核なしデータの読み込み
        file_no = file_info(idx_no_nucleus(j)).filename;
        [OOP_no, grid_no] = calculate_OOP_and_fiber_orientation_grid(file_no);

        % 核ありデータの読み込み
        file_with = file_info(idx_with_nucleus(j)).filename;
        [OOP_with, grid_with] = calculate_OOP_and_fiber_orientation_grid(file_with);

        % COOPの計算
        COOP = calculate_COOP(OOP_no, OOP_with, grid_no, grid_with);
        % 核の有無を結果に記録
        COOP_results = [COOP_results; struct('final_value', fv, 'nucleus_obstruction', 0, 'COOP', COOP)];
    end
end

% 結果をデータテーブルに変換
COOP_table = struct2table(COOP_results);

% データをソート
[sorted_final_values, sort_idx] = sort(COOP_table.final_value);
sorted_COOP = COOP_table.COOP(sort_idx);
sorted_nucleus_obstruction = COOP_table.nucleus_obstruction(sort_idx);

% グループ変数をカテゴリカルに変換
FinalValue = categorical(sorted_final_values);

% ANOVAによる統計検定とTukey HSD（NucleusObstructionを除外）
group = FinalValue;
[p, tbl, stats] = anovan(sorted_COOP, sorted_final_values, 'varnames', {'FinalValue'}, 'display', 'off');

% Tukey HSD (多重比較)
[c, m, h, gnames] = multcompare(stats, 'CType', 'hsd', 'Display', 'off');


% ---- グラフをプロット ----
% 'No Obstruction' と 'With Obstruction' のデータを統合し、平均と標準偏差を計算
unique_final_values = categories(FinalValue);
unique_final_values_num = str2double(unique_final_values);

COOP_mean = zeros(length(unique_final_values),1);
COOP_sd = zeros(length(unique_final_values),1);

for i = 1:length(unique_final_values)
    fv = unique_final_values{i};
    idx = FinalValue == fv;
    COOP_mean(i) = mean(sorted_COOP(idx));
    COOP_sd(i) = std(sorted_COOP(idx));
end

% バープロットの作成
figure;
bar(unique_final_values_num, COOP_mean, 'FaceColor', 'b', 'BarWidth', 0.5);
hold on;
errorbar(unique_final_values_num, COOP_mean, COOP_sd, '.', 'Color', 'k');
xlabel('Final Value (Stretch Condition)');
ylabel('Mean COOP');
ylim([0, 1]);
title('COOP for Different Final Values');
grid on;

% 有意差を表示（p < 0.05 のもののみ）
% 有意差表示用の y 軸の開始位置と増分を設定
y_base = max(COOP_mean + COOP_sd) + 0.1; % ラインの開始位置を設定
y_increment = 0.05; % ラインの高さの増分

% 現在のラインの数を記録するカウンターを初期化
line_counter = 0;

for i = 1:length(unique_final_values)
    for j = i+1:length(unique_final_values)
        % グループ番号を取得
        group1 = find(strcmp(gnames, ['FinalValue=' unique_final_values{i}]));
        group2 = find(strcmp(gnames, ['FinalValue=' unique_final_values{j}]));
        
        % 対比較を探す
        comp_idx = ( (c(:,1) == group1 & c(:,2) == group2) | (c(:,1) == group2 & c(:,2) == group1) );
        if any(comp_idx)
            p_value = c(comp_idx, 6);
            if p_value < 0.05
                % 有意差をプロット
                x1 = unique_final_values_num(i) + 0.1;
                x2 = unique_final_values_num(j) - 0.1;
                y = y_base + line_counter * y_increment;
                plot([x1, x2], [y, y], '-k', 'LineWidth', 1.5);
                text(mean([x1, x2]), y + 0.02, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
                % カウンターを増加
                line_counter = line_counter + 1;
            end
        end
    end
end


hold off;

% ---- 結果を表示 ----
disp('ANOVA results for COOP:');
disp(tbl);
disp('Tukey HSD multiple comparison results:');
disp(array2table(c, 'VariableNames', {'Group1','Group2','LowerLimit','MeanDiff','UpperLimit','pValue'}));


%% 関数定義

% 関数: OOPとグリッドごとの繊維方向を計算
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

% 関数: COOPを計算
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

% 関数: グリッドごとの繊維方向ベクトルを収集
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

% 関数: OOPを計算
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