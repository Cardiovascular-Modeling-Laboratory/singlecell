% The purpose of this function is to stretch the cell geometry horizontally for the fixed amount based on the previous cell geometry
% For now, only square cell is considered
%
function [mat_r,Npts_t,drx,dry,dA,dr_dist_squared, dist_to_line_sq, shape_name,dist_pair,drx_norm,dry_norm,Concave_ind,outline,nuc_x,nuc_y,nuc_cx,nuc_cy,choice, outside_segs, inside_segs,bdry_mat,bdry_pts,out_ind,Num_points,A,nuc_radius]=cell_geometry_strecth_v1(A_initial, stretch_factor, mat_r_initial, Npts_t,drx,dry,dA,dr_dist_squared, dist_to_line_sq, shape_name,dist_pair,drx_norm,dry_norm,Concave_ind,outline_initial,nuc_x_initial,nuc_y,nuc_cx_initial,nuc_cy,choice, outside_segs, inside_segs,bdry_mat,bdry_pts,out_ind,Num_points,A,nuc_radius, poisson_ratio)


% この関数の目的は、既存のセルジオメトリを基に水平方向に一定量だけセルを伸ばすことです。
% 現在のところ、正方形のセルのみを考慮します。
% 伸縮係数を設定（例：1.5倍に伸ばす）
% Stretch factor (integer)

% Poisson ratioを設定
% mat_r の x 座標を伸ばす
mat_r_initial(:,1) = mat_r_initial(:,1) * stretch_factor;
mat_r = mat_r_initial;
% mat_r　の y 座標をPoission ratioに基づいて伸ばす
mat_r_initial(:,2) = mat_r_initial(:,2) * (1 + poisson_ratio * (stretch_factor - 1));

% セルの面積 A を更新
A = A_initial * stretch_factor^2;

% outline の x 座標を伸ばす
outline_initial(1,:) = outline_initial(1,:) * stretch_factor;
% outline の y 座標をPoission ratioに基づいて伸ばす
outline_initial(2,:) = outline_initial(2,:) * (1 + poisson_ratio * (stretch_factor - 1));
outline = outline_initial;


% 核の位置と形状を更新
% initial使わないと、指数的増加するぞ
nuc_cx = nuc_cx_initial * stretch_factor;
nuc_x = nuc_x_initial * stretch_factor;

% x 座標のユニークな値を取得
mtx = unique(mat_r(:,1));
mty = unique(mat_r(:,2));

% ポイント間の最小距離を計算
xd = min(abs(mtx(1:end-1) - mtx(2:end))) / 2;
yd = min(abs(mty(1:end-1) - mty(2:end))) / 2;

% r_col_x, r_row_x, r_col_y, r_row_y を再計算
r_col_x = repmat(mat_r(:,1), 1, Npts_t);
r_row_x = r_col_x';
r_col_y = repmat(mat_r(:,2), 1, Npts_t);
r_row_y = r_col_y';

% drx, dry を再計算
drx = r_row_x - r_col_x;
dry = r_row_y - r_col_y;

% 距離の二乗を再計算
dr_dist_squared = drx.^2 + dry.^2;
dist_pair = sqrt(dr_dist_squared);

% 正規化された drx, dry を再計算
drx_norm = drx ./ dist_pair;
drx_norm(isnan(drx_norm)) = 0;
dry_norm = dry ./ dist_pair;
dry_norm(isnan(dry_norm)) = 0;

% dist_to_line_sq を再計算
N = Npts_t;
drx_trip_constr3 = repmat(drx, [1, 1, N]);
dry_trip_constr3 = repmat(dry, [1, 1, N]);
drx_trip_constr = permute(drx_trip_constr3, [3, 2, 1]);
dry_trip_constr = permute(dry_trip_constr3, [3, 2, 1]);
dr_dist_squaredr = permute(repmat(dr_dist_squared, [1, 1, N]), [3, 2, 1]);
dot_product_trip = drx_trip_constr .* drx_trip_constr3 + dry_trip_constr .* dry_trip_constr3;
to_line_x = drx_trip_constr3 .* dr_dist_squaredr - drx_trip_constr .* dot_product_trip;
to_line_y = dry_trip_constr3 .* dr_dist_squaredr - dry_trip_constr .* dot_product_trip;
dist_to_line_sq_temp = (to_line_x.^2 + to_line_y.^2) ./ (dr_dist_squaredr.^2);

% NaN をゼロに置き換える
temp_zero = zeros(size(dist_to_line_sq_temp));
dist_to_line_sq = dist_to_line_sq_temp;
dist_to_line_sq(find(isnan(dist_to_line_sq_temp))) = temp_zero(find(isnan(dist_to_line_sq_temp)));
clear dist_to_line_sq_temp temp_zero
%normalized drx and dry and get rid of NaN
temp_zero = zeros(size(drx));
drx_norm = drx./dist_pair;
drx_norm(find(isnan(drx_norm))) = temp_zero(find(isnan(drx_norm)));
temp_zero = zeros(size(dry));
dry_norm = dry./dist_pair;
dry_norm(find(isnan(dry_norm))) = temp_zero(find(isnan(dry_norm)));

% Concave_ind を設定（正方形なので 1）
Concave_ind = 1;

% 形状名を更新
shape_name = 'square';

% 単位面積 dA を再計算
if length(Concave_ind) ~= 1
    dA = A / sum(Concave_ind);
else
    dA = A / Npts_t;
end

% 境界セグメントを再計算
inside_segs = ones(size(drx));
outside_segs = zeros(size(drx));
bdry_mat = zeros(size(drx));

% 境界点を再計算
bdry_pts = [];
geo = outline';
t = linspace(0, 1);

for j=1:size(geo,1)-1
    P1=geo(j,:);
    P2=geo(j+1,:);
    seg_x=P1(1).*t+P2(1).*(1-t);
    seg_y=P1(2).*t+P2(2).*(1-t);
    Seg_x1=[seg_x'-xd seg_y'];
    Seg_x2=[seg_x'+xd seg_y'];
    Seg_y1=[seg_x' seg_y'-yd];
    Seg_y2=[seg_x' seg_y'+yd];
    Box=[Seg_x1;Seg_x2; Seg_y1; Seg_y2];
    k=convhull(Box(:,1),Box(:,2));
    Region=Box(k,:);
    in=inpoly(mat_r,Region);
    bdry_pt_idx=find(in==1);
    bdry_pts=[bdry_pts; mat_r(bdry_pt_idx,:)];
    
    meh=[nchoosek(bdry_pt_idx,2);fliplr(nchoosek(bdry_pt_idx,2))];
    bdry_mat(sub2ind(size(bdry_mat),meh(:,1),meh(:,2)))=1;
end
% [out_ind,bdry_pts] = boundary_determ_v2(outline, mat_r,choice,bdry_pts);

% Convex Hull を使用して境界点を整理
out_ind = boundary(mat_r(:,1), mat_r(:,2));
bdry_pts = mat_r(out_ind, :);

% 必要に応じてプロット（オプション）
% figure;
% plot(mat_r(:,1), mat_r(:,2), '*r');
% hold on;
% plot(bdry_pts(:,1), bdry_pts(:,2), 'b-', 'LineWidth', 2); % 境界線を太くして見やすくする
% plot(nuc_x, nuc_y, 'g', 'LineWidth', 1.5); % 核の位置も強調
% grid on;
% hold off;

% title(['Cell Geometry with Stretch Factor = ', num2str(stretch_factor)]);
% xlabel('X Position');
% ylabel('Y Position');

% % drx, dry, drx_norm, dry_normのplot
% figure;
% subplot(2,2,1);
% plot(drx(:), dry(:), '*r');
% title('drx vs dry');
% xlabel('drx');
% ylabel('dry');
% grid on;
% subplot(2,2,2);
% plot(drx_norm(:), dry_norm(:), '*r');
% title('drx_{norm} vs dry_{norm}');
% xlabel('drx_{norm}');
% ylabel('dry_{norm}');
% grid on;
% subplot(2,2,3);
% plot(drx(:), drx_norm(:), '*r');
% title('drx vs drx_{norm}');
% xlabel('drx');
% ylabel('drx_{norm}');
% grid on;
% subplot(2,2,4);
% plot(dry(:), dry_norm(:), '*r');
% title('dry vs dry_{norm}');
% xlabel('dry');
% ylabel('dry_{norm}');
% grid on;

end
