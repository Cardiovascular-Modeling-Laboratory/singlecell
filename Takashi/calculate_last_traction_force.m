function max_force = calculate_last_traction_force(F_store)
    % 最大牽引力を計算する関数
    % F_store: 力のデータ (Nx2 の行列、[F_x, F_y])
    % max_force: 最大牽引力

    % 各時点での牽引力の大きさを計算
    forces = sqrt(F_store(:, :, 1).^2 + F_store(:, :, 2).^2); % 合力
    % max_force = max(forces(:)); % 最大値を取得
    % use last time step force
    % 絶対値の最大値を取得
    max_force = max(forces(end, :));
end
