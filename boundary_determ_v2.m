% this code determines the boundary points of the cell geometry
function [out_ind,bdry_pts]=boundary_determ_v2(outline, mat_r,choice,bdry_pts)

if choice==3
    dist_mat=sqrt(mat_r(:,1).^2+mat_r(:,2).^2); % distance between each lattice point ant the center of the circle
    % since we want to isolate the boundary points of the lattice,
    % allow for a little wiggle room (for roundoff error) when
    % determining which distances are "approximately equal" to the
    % circle radius
    out_ind_temp=dist_mat>=0.99.*outline(3);
    out_ind=find(out_ind_temp==1);
elseif choice==7
    outline_in=0.99.*outline; % inner oval (oval just slightly inside of the outline)
    outline_out=1.01.*outline; % outter oval (oval just slightly outside of the outline)
    out_ind_temp=zeros(size(mat_r,1),1);
    dist_mat=sqrt(mat_r(:,1).^2+mat_r(:,2).^2); % distance of each lattice point to center of ellipse
    % for each lattice point, construct a ray pointing from the oval's
    % center to a point outside the outter oval. then, determine where
    % this ray interests with the inner and outer ovals. compute the
    % distances of these 2 points. if the lattice point distance is
    % between these 2 computed distances, then it is on the boundary:
    for j=1:size(mat_r,1)
        if sign(mat_r(j,1))==0 % determine if the lattice point is along the vertical axis
            % construct ray
            lx=zeros(1,100);
            ly=linspace(0,3*(0.5+sign(mat_r(j,2)))*max(mat_r(:,2)));
        else
            % determine which quadrant the point is in
            X=mat_r(j,1)>=0;
            Y=mat_r(j,2)>=0;
            quad_number=3+X-Y-2*X*Y;
            % depending on the quadrant, impose the proper sign
            switch quad_number
                case 1 
                    c=1;
                case 2 
                    c=-1;
                case 3
                    c=-1;
                case 4
                    c=1;
            end
            % construct ray
            xmax=2*c*max(mat_r(:,1));
            lx=linspace(0,xmax);
            ly=linspace(0,(mat_r(j,2)/mat_r(j,1))*xmax);
        end
        [x_in, y_in]=intersections(lx,ly,outline_in(1,:),outline_in(2,:)); % determines where ray interestcs the inner ellipse
        [x_out, y_out]=intersections(lx,ly,outline_out(1,:),outline_out(2,:)); % determines where ray interestcs the outer ellipse
        dist_in=sqrt(x_in.^2+y_in.^2); % distance between inner ntersection point and the center
        dist_out=sqrt(x_out.^2+y_out.^2); % distance between outer ntersection point and the center
        out_ind_temp(j)=((dist_mat(j)>= dist_in)&(dist_mat(j)<=dist_out)); % determine if the lattice point distance is within the distance range
    end
    out_ind=find(out_ind_temp==1);
else
    out_ind=find(ismember(mat_r,bdry_pts,'rows')==1);
end

bdry_pts=mat_r(out_ind,:);

end % end function

