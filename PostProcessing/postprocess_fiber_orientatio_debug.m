% fiber_multiple_orientation_visualization_with_curves.m
% This script loads a simulation result file, extracts the coordinates of three fibers,
% calculates their angles, and plots them as curves to visualize if the angle calculation is correct.

% 読み込むシミュレーション結果ファイル名（例として最初のファイルを選択）
filePattern = fullfile(pwd, 'simulation_results_*.mat');
resultFiles = dir(filePattern);
filename = fullfile(resultFiles(1).folder, resultFiles(1).name);  % 1番目のファイルを使用
load(filename, 'net_res', 'nuc_cx_store', 'nuc_cy_store', 'nuc_radius');

% シミュレーション条件の設定
j2 = 1; % nucleus placement index
m = 1; % 1つ目のシミュレーション
net1 = net_res{j2, m};
net = net1{end}; % 最終結果のみを使用

% 可視化するための3つのfiberを選択
fiber_indices = [1, 2, 3]; % 3つのfiberのインデックス
fiber_segment_index = 1; % 各fiberの1つ目のsegmentを使用

% カラーマップの設定
colors = lines(3);

figure;
hold on;
title('Multiple Fiber Orientation Visualization with Curves', 'FontSize', 20);
xlabel('X Coordinate', 'FontSize', 20);
ylabel('Y Coordinate', 'FontSize', 20);
axis equal;
grid on;

% 各fiberをプロットし、角度を表示
for i = 1:length(fiber_indices)
    fiber_index = fiber_indices(i);
    net_temp = net{fiber_index};
    
    % fiberの全ての点を取得してプロット
    x_coords = net_temp(1, :, fiber_segment_index);
    y_coords = net_temp(2, :, fiber_segment_index);
    plot(x_coords, y_coords, 'Color', colors(i, :), 'LineWidth', 1.5); % 曲線としてfiberをプロット
    
    % fiberの始点と終点を使って角度を計算
    x_start = x_coords(1);
    y_start = y_coords(1);
    x_end = x_coords(end);
    y_end = y_coords(end);
    dx = x_end - x_start;
    dy = y_end - y_start;
    angle = atan2d(dy, dx); % 度数法で角度を取得
    
    % 90度を超える場合は 180 - 角度 に変換（0-90度に収める）
    if angle > 90
        angle = 180 - angle;
    end
    
    % 方向を示す矢印のプロット
    quiver(x_start, y_start, dx, dy, 0, 'Color', colors(i, :), 'MaxHeadSize', 1, 'LineWidth', 1.5);
    
    % 角度を表示
    mid_x = (x_start + x_end) / 2;
    mid_y = (y_start + y_end) / 2;
    text(mid_x, mid_y, sprintf('%.1f°', angle), 'FontSize', 20, ...
         'HorizontalAlignment', 'center', 'Color', colors(i, :));
end

hold off;
