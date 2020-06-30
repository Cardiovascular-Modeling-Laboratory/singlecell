% input vector v and determine the 5 vectors that will sum up to v. Here,
% alpha(k) is the magnitude of each decomposition vector w_k 
% input the number of fibers in the decomposition (N_F), the maximum
% allowed number of fibers in the decomposition (N_F_max), the lattice
% point of interest from mat_r (pt), the cell outline (outline), and the
% lattice (mat_r)
% See notes for more explicit derivation/explanation
% Output a matrix of vectors W and an vale W_n specifying the number
% of vectors used in the decomposition 

function [W, W_n] = v_splitting3(v,N_F,N_F_max,pt,outline,Lat,mat_r,choice,out_ind,outside_segs,lat_sat_idx)
% the construction requires v=[v_x v_y] to be in quadrant 1: v_x,v_y>0. So
% given v, use the transformation matrix T to create a new vector
% satisfying this condition
v_normalized=v./norm(v); % normalized v
v_x=v_normalized(1);
v_y=v_normalized(2);

T=[sign(v_x) 0; 0 sign(v_y)]; % transformation matrix
v_Q1=v_normalized*T; % the transformed and normalized v

%% given the lattice point and force vector, determine the relevant angles for the decomposition
% identify the lattice points which lie on the cell boundary
%out_ind=boundary_determ(outline, mat_r,choice);
bdry_lat=Lat(out_ind,:);

% determine if the point of interest lies on the boundary or the interior
% of the cell
on_cell=ismember(pt,bdry_lat,'rows'); % on_cell=1 if point is a boundary point, 0 if point is on interior

% determine the points on the concave portion of the cell that are
% saturated
lat_sat_idx_concave=intersect(find(sum(outside_segs)>0),lat_sat_idx);
bdry_lat_concave=Lat(lat_sat_idx_concave,:);
on_cell_concave=ismember(pt,bdry_lat_concave,'rows');

switch on_cell_concave
    case 1 % the case where the point of interest is a boundary point on the concave region of the cell
        % find two nearest boundary points to p1_bdry
        p1_bdry_idx=find(ismember(bdry_lat,pt,'rows'));
        bdry_lat_temp=bdry_lat;
        bdry_lat_temp(p1_bdry_idx,:)=[];
        q_pts=zeros(2,2);
        for j=1:2
            p1_bdry_rep=repmat(pt,size(bdry_lat_temp,1),1);
            p1_bdry_dist=sqrt((bdry_lat_temp(:,1)-p1_bdry_rep(:,1)).^2+(bdry_lat_temp(:,2)-p1_bdry_rep(:,2)).^2);
            [~,min_idx]=min(p1_bdry_dist);
            q_i=bdry_lat_temp(min_idx(1),:);
            q_pts(j,:)=q_i;

            % now that the relevant point has been identified, need to identify new
            % points so remove the identified point from the point list
            qi_idx=find(ismember(bdry_lat_temp,q_i,'rows'));
            bdry_lat_temp(qi_idx,:)=[];
        end
        % construct vectors that'll be used to aid in determining theta_max
        v1=q_pts(1,:)-pt;
        v2=q_pts(2,:)-pt;
        
        theta_F1=acos(dot(v1*T,v_Q1)/sqrt(sum((v1*T).^2)*sum(v_Q1.^2)));
        theta_F2=acos(dot(v2*T,v_Q1)/sqrt(sum((v2*T).^2)*sum(v_Q1.^2)));

        % determine the angle between v1 and v2
        theta_v1_v2=acos(dot(v1,v2)/sqrt(sum(v1.^2)*sum(v2.^2)));
        
        % if this angle is larger than, say, 150 degrees then take
        % theta_max to be the smaller of the 2 potential angles. otherwise,
        % take the larger of the 2 potential angles
        V12=[v1; v2];
        if theta_v1_v2>150*pi/180
            theta_max_temp=min(theta_F1,theta_F2);
        else
            theta_max_temp=max(theta_F1,theta_F2);
        end
        temp_idx=find([theta_F1,theta_F2]==theta_max_temp);
        v_min=V12(temp_idx,:);
        
        % now bisect the angle formed by v and v_min
        v_bisect=norm(v).*v_min+norm(v_min).*v;
        
        % the 3 vectors to be considered for the remaining portion of the
        % code are now the original vector v, the bisector v_bisect, and
        % the vector v_min. Relabel them
        v1=v_min;
        v2=v;
        v=v_bisect;
        v_normalized=v./norm(v); % normalized v, redefined
        v_x=v_normalized(1);
        v_y=v_normalized(2);
        T=[sign(v_x) 0; 0 sign(v_y)]; % transformation matrix
        v_Q1=v_normalized*T; % the transformed and normalized v

    otherwise
        switch on_cell
            case 1 % the case when the point of interest is a boundary point but not on the concave region of the cell
                % find two nearest boundary points to p1_bdry
                p1_bdry_idx=find(ismember(bdry_lat,pt,'rows'));
                bdry_lat_temp=bdry_lat;
                bdry_lat_temp(p1_bdry_idx,:)=[];
                q_pts=zeros(2,2);
                for j=1:2
                    p1_bdry_rep=repmat(pt,size(bdry_lat_temp,1),1);
                    p1_bdry_dist=sqrt((bdry_lat_temp(:,1)-p1_bdry_rep(:,1)).^2+(bdry_lat_temp(:,2)-p1_bdry_rep(:,2)).^2);
                    [~,min_idx]=min(p1_bdry_dist);
                    q_i=bdry_lat_temp(min_idx(1),:);
                    q_pts(j,:)=q_i;

                    % now that the relevant point has been identified, need to identify new
                    % points so remove the identified point from the point list
                    qi_idx=find(ismember(bdry_lat_temp,q_i,'rows'));
                    bdry_lat_temp(qi_idx,:)=[];
                end
                % construct the vectors that will be used to determine the max
                % angle. 
                v1=q_pts(1,:)-pt;
                v2=q_pts(2,:)-pt;
            otherwise % the case when the point of interest is an interior point
                % construct the vectors that will be used to determine the max
                % angle. since this is an interior point, we'll construct vectors
                % which are orthogonal to F_est by rotating F_est by plus/minus 90
                % degrees
                theta_plus=pi/2;
                theta_minus=-pi/2;
                R_plus=[cos(theta_plus), -sin(theta_plus); sin(theta_plus), cos(theta_plus)];
                R_minus=[cos(theta_minus), -sin(theta_minus); sin(theta_minus), cos(theta_minus)];

                v1=(R_plus*v_normalized')';
                v2=(R_minus*v_normalized')';
        end
end

% determine the angle between vectors v1 & v_Q1 and between v2 & v_Q1. Note
% that we need to multiply v1 and v2 by the transformation matrix T since
% v_Q1 is the transformed version of v
theta_F1=acos(dot(v1*T,v_Q1)/sqrt(sum((v1*T).^2)*sum(v_Q1.^2)));
theta_F2=acos(dot(v2*T,v_Q1)/sqrt(sum((v2*T).^2)*sum(v_Q1.^2)));


% record theta_F and theta_max. theta_F is the angle v_Q1 makes with the
% positive x-axis and theta_max is the maximum allowable angle in the
% decomposition. Enforce the restriction that the angle between w1 and w5
% can be no more than 170 degrees, i.e., theta_max can be no more than 85
% degrees
theta_F=atan(v_Q1(2)/v_Q1(1));
theta_max=min(min(theta_F1,theta_F2),85*(pi/180));

%% perform the decomposition based on how many vectors will exist in the decomposition

alpha=zeros(1,N_F_max);

k=1:N_F_max;
th=ones(1,5);
th([2 4])=0.5;
theta=theta_F+sign(3-k).*th.*theta_max;

w_x=cos(theta);
w_y=sin(theta);

% complete explanations of where these expressions come from should be in
% the notes
switch N_F
    case 1
        alpha(3)=1;
    case 3
        alpha(2)=1/(2*cos(theta_max/2)+3);
        alpha(3)=3*alpha(2);
        alpha(4)=alpha(2);
    case 5
        alpha(1)=1/(2*cos(theta_max)+4*cos(theta_max/2)+4);
        alpha(2)=2*alpha(1);
        alpha(3)=2*alpha(2);
        alpha(4)=alpha(2);
        alpha(5)=alpha(1);
    otherwise
        disp(['Number of decomposition vectors is not acceptable!:', num2str(N_F)])
        return
end



% construct w as a matrix
W_x=alpha.*w_x;
W_y=alpha.*w_y;
W_Q1=[W_x; W_y]'; % the matrix of decomposition vectors. 

% scale each w to be on par with original v
W_temp=W_Q1.*norm(v);

% we now have the vector decomposition with correct scaling. All that's
% left to do is transform it so that it lines up with the original position
% of v. So multiply each decomposed vector by the transformation matrix T
% (multiply by T since T is it's own inverse)
W=zeros(N_F_max,2);
for j=1:N_F_max
    W(j,:)=W_temp(j,:)*T;
end

% right now, W may contain rows of 0s since some alpha's may be 0 (when the
% number of decomposition vectors is 1 or 3). in this case, we don't care
% about the other rows. so only hold on to the rows that are relevant
switch N_F
    case 1
        W=W(3,:);
    case 3
        W=W(2:4,:);
end
W_n=size(W,1);




