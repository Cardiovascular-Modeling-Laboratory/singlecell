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
