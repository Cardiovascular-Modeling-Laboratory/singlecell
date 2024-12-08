% plot_fiber_orientation_comparison_subplots.m
% This script loads all simulation results files from the current directory,
% calculates the fiber network orientation angles (with adjustment to keep angles within [0, 90] degrees),
% and creates a comparison plot with each max_stretch_factor shown in a separate subplot.

% ファイル名のパターンに基づいてシミュレーション結果を集める
filePattern = fullfile(pwd, 'simulation_results_*.mat');
resultFiles = dir(filePattern);

% max_stretch_factor 別に角度データを保存するための構造体を初期化
angleData = struct();

% 各シミュレーション結果ファイルに対してループ
for i = 1:length(resultFiles)
    % ファイル名から max_stretch_factor を抽出
    filename = resultFiles(i).name;
    max_stretch_factor = sscanf(filename, 'simulation_results_%f');
    
    % max_stretch_factor を文字列化してフィールド名に変換（不正文字を削除）
    factorStr = matlab.lang.makeValidName(sprintf('factor_%.1f', max_stretch_factor));

    % ファイルの読み込み
    fullFilePath = fullfile(resultFiles(i).folder, filename);
    load(fullFilePath, 'net_res', 'nuc_cx_store', 'nuc_cy_store', 'nuc_radius');
    
    % シミュレーション条件の設定
    j2 = 1; % nucleus placement index
    nsim = size(net_res, 2);
    
    % 角度配列を初期化（この max_stretch_factor 用）
    angles = [];

    % 各シミュレーションに対してループ
    for m = 1:nsim
        % fiber network のデータを取得
        net1 = net_res{j2, m};
        net = net1{end}; % 最終結果のみを使用

        % 繊維の角度を計算
        for j = 1:size(net, 1)
            net_temp = net{j};
            for k = 1:size(net_temp, 3)
                % 繊維のx座標とy座標を取得
                net_x_temp = net_temp(1, :, k);
                net_y_temp = net_temp(2, :, k);

                % 始点と終点の差を計算
                dx = net_x_temp(end) - net_x_temp(1);
                dy = net_y_temp(end) - net_y_temp(1);

                % 角度を計算
                angle = atan2d(dy, dx); % 度数法で角度を取得
                
                % 90度を超える場合は 180 - 角度 に変換
                if angle > 90
                    angle = 180 - angle;
                end
                
                % リストに追加
                angles = [angles; angle];
            end
        end
    end

    % 角度データを構造体に保存
    if isfield(angleData, factorStr)
        angleData.(factorStr) = [angleData.(factorStr); angles];
    else
        angleData.(factorStr) = angles;
    end
end

% すべての角度データから分布の最大 y 値を取得して ylim を統一
fields = fieldnames(angleData);
numFactors = length(fields);
yLimits = 0; % 初期値
for i = 1:numFactors
    angles = angleData.(fields{i});
    h = histogram(angles, 'BinEdges', 0:5:90, 'Normalization', 'probability', 'Visible', 'off');
    yLimits = max(yLimits, max(h.Values));
    
    % ヒストグラムの頻度の合計を確認
    freq_sum = sum(h.Values);
    if abs(freq_sum - 1.0) > 1e-6
        warning('Histogram frequency sum for %s is not 1.0 (sum = %.4f)', fields{i}, freq_sum);
    end
end

% 各 max_stretch_factor の角度分布をサブプロットで表示
figure;
for i = 1:numFactors
    subplot(ceil(sqrt(numFactors)), ceil(sqrt(numFactors)), i);
    
    % 各角度データをプロット
    angles = angleData.(fields{i});
    histogram(angles, 'BinEdges', 0:5:90, 'Normalization', 'probability');
    xlabel('Angle (degrees)');
    ylabel('Frequency');
    title(['max_stretch_factor = ', strrep(fields{i}, 'factor_', '')], 'Interpreter', 'none'); % LaTeX オフ
    ylim([0 yLimits]);
    grid on;
end

% 全体の図として保存
saveas(gcf, 'FiberOrientationComparison_Subplots.png');
