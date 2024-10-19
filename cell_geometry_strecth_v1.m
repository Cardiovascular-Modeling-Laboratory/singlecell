% The purpose of this function is to stretch the cell geometry horizontally for the fixed amount based on the previous cell geometry
% For now, only square cell is considered
%
function [mat_r, Npts_t, drx, dry, dA, dr_dist_squared, dist_to_line_sq, shape_name, dist_pair, drx_norm, dry_norm, Concave_ind, outline, nuc_x, nuc_y, nuc_cx, nuc_cy, choice, outside_segs, inside_segs, bdry_mat, bdry_pts, out_ind, Num_points, A, nuc_radius] = cell_geometry_stretch_v1(A_initial, stretch_factor, mat_r_initial, Npts_t, drx, dry, dA, dr_dist_squared, dist_to_line_sq, shape_name, dist_pair, drx_norm, dry_norm, Concave_ind, outline_initial, nuc_x_initial, nuc_y_initial, nuc_cx_initial, nuc_cy_initial, choice, outside_segs, inside_segs, bdry_mat, bdry_pts, out_ind, Num_points, A, nuc_radius, poisson_ratio)

    % Stretch in x-direction
    mat_r_initial(:,1) = mat_r_initial(:,1) * stretch_factor;

    % Shrink in y-direction based on Poisson's ratio
    mat_r_initial(:,2) = mat_r_initial(:,2) / (1 + poisson_ratio * (stretch_factor - 1));

    mat_r = mat_r_initial;

    % Update area based on stretch in x and shrink in y
    A = A_initial * stretch_factor * (1 / (1 + poisson_ratio * (stretch_factor - 1)));

    % Stretch the outline in x-direction
    outline_initial(1,:) = outline_initial(1,:) * stretch_factor;

    % Shrink the outline in y-direction based on Poisson's ratio
    outline_initial(2,:) = outline_initial(2,:) / (1 + poisson_ratio * (stretch_factor - 1));

    outline = outline_initial;

    % Update nuclear position for both x and y directions
    % Stretch the x position of the nucleus
    nuc_cx = nuc_cx_initial * stretch_factor;
    nuc_x = nuc_x_initial * stretch_factor;

    % Shrink the y position of the nucleus based on Poisson's ratio
    nuc_cy = nuc_cy_initial / (1 + poisson_ratio * (stretch_factor - 1));
    nuc_y = nuc_y_initial / (1 + poisson_ratio * (stretch_factor - 1));


    % Unique x and y points
    mtx = unique(mat_r(:,1));
    mty = unique(mat_r(:,2));

    % Calculate minimum distance between points
    xd = min(abs(mtx(1:end-1) - mtx(2:end))) / 2;
    yd = min(abs(mty(1:end-1) - mty(2:end))) / 2;

    % Recalculate drx, dry
    r_col_x = repmat(mat_r(:,1), 1, Npts_t);
    r_row_x = r_col_x';
    r_col_y = repmat(mat_r(:,2), 1, Npts_t);
    r_row_y = r_col_y';

    drx = r_row_x - r_col_x;
    dry = r_row_y - r_col_y;

    % Recalculate squared distances
    dr_dist_squared = drx.^2 + dry.^2;
    dist_pair = sqrt(dr_dist_squared);

    % Normalize drx and dry
    drx_norm = drx ./ dist_pair;
    dry_norm = dry ./ dist_pair;

    drx_norm(isnan(drx_norm)) = 0;
    dry_norm(isnan(dry_norm)) = 0;

    % Calculate dist_to_line_sq (if needed, keep as before)

    % Set Concave_ind (for square, it's 1)
    Concave_ind = 1;

    % Update shape name
    shape_name = 'square';

    % Recalculate unit area dA
    dA = A / Npts_t;

    % Recalculate boundary points (same as original)
    bdry_pts = [];
    geo = outline';
    t = linspace(0, 1);

    for j = 1:size(geo,1)-1
        P1 = geo(j,:);
        P2 = geo(j+1,:);
        seg_x = P1(1).*t + P2(1).*(1-t);
        seg_y = P1(2).*t + P2(2).*(1-t);
        Seg_x1 = [seg_x'-xd seg_y'];
        Seg_x2 = [seg_x'+xd seg_y'];
        Seg_y1 = [seg_x' seg_y'-yd];
        Seg_y2 = [seg_x' seg_y'+yd];
        Box = [Seg_x1; Seg_x2; Seg_y1; Seg_y2];
        k = convhull(Box(:,1), Box(:,2));
        Region = Box(k,:);
        in = inpoly(mat_r, Region);
        bdry_pt_idx = find(in == 1);
        bdry_pts = [bdry_pts; mat_r(bdry_pt_idx,:)];

        meh = [nchoosek(bdry_pt_idx,2); fliplr(nchoosek(bdry_pt_idx,2))];
        bdry_mat(sub2ind(size(bdry_mat), meh(:,1), meh(:,2))) = 1;
    end

    % Update boundary points and indices
    out_ind = boundary(mat_r(:,1), mat_r(:,2));
    bdry_pts = mat_r(out_ind, :);

end
