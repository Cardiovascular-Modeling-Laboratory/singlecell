% run_simulations_with_OOP_calculation.m
% 複数の max_stretch_factor と fold_time でシミュレーションを実行し、
% Orientation Order Parameter を計算します。
% エラーが発生した場合、エラーログを保存し、次のシミュレーションを継続します。

% max_stretch_factors と fold_time の設定
% max_stretch_factors = [3.5, 4.0, 5.0]; % 使用する max_stretch_factor のリスト
% fold_times = [(72-1)*3600, (72-3)*3600, (72-6)*3600, (72-12)*3600, (72-24)*3600]; % 使用する fold_time のリスト
% fold_times = [(72-1)*3600]; % 使用する fold_time のリスト
max_stretch_factors = [1.0, 2.0, 3.0, 4.0, 5.0]; % 使用する max_stretch_factor のリスト
nuc_rel_vec=[0, 1];
n_trial = 5; % 各条件での試行回数

% エラーログ用の保存ディレクトリ
error_log_dir = fullfile(pwd, 'error_logs');
if ~exist(error_log_dir, 'dir')
    mkdir(error_log_dir);
end

% 各 max_stretch_factor でループを回してシミュレーションを実行
for i = 1:length(max_stretch_factors)
    max_stretch_factor = max_stretch_factors(i);
    for j = 1:length(nuc_rel_vec)
        nuc_rel = nuc_rel_vec(j);

        for k = 1:n_trial
            try
                % single_cell_units_linked_v3 を呼び出し、結果を保存
                disp(['Running simulation with max_stretch_factor = ', num2str(max_stretch_factor), ...
                    ' and nuc_rel = ', num2str(nuc_rel)]);
                single_cell_units_linked_v3_NucPlacement(max_stretch_factor, [nuc_rel], k);

                % Orientation Order Parameter を計算
                % calculate_orientation_order_parameter;
                
            catch ME
                % エラー情報を保存
                errorTime = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
                errorMessage = ME.message;
                errorDetails = struct('Time', errorTime, 'Message', errorMessage, 'Stack', ME.stack);
                error_filename = fullfile(error_log_dir, ['error_', errorTime, '_factor_', ...
                                num2str(max_stretch_factor), '_fold_', num2str(nuc_rel), '.mat']);
                save(error_filename, 'errorDetails');
                disp(['Error occurred and saved to ', error_filename]);
                disp(['Error message: ', errorMessage]);
                disp(['Error stack: ']);
                ME.stack
            end
        end
    end
    close all;
end
