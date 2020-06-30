function E_BD=E_BD_curve(rtx_test,rty_test,boundary_pts,t,dist_pair,P0,P1,P2,P3,P4,dA_FA_n)
k_stretch=0.3; % ~5 k_B*T/um^2 (fitted parameter)
A_FA=dA_FA_n; % focal adhesion area for nascent focal adhesion in um^2 (this should be the same value used in the force equation code but with converted units since all values in the fiber code are in um scale)
% since the curve is constructed at the micrometer scale (k_B*T is factored
% out in the energy equation so it's not here)

%% determine curve indices which travel outside the boundary
test_curve = [rtx_test' rty_test'];
in_bdry = inpoly(test_curve,boundary_pts);
% due to possible numerical error, we'll set the starting and ending
% points of in_bdry to 0. This way, the starting and ending points
% of the curve must be considered
in_bdry([1, 100])=[];
s=find(in_bdry==0); % indices along the curve where the curve passes outside the boundary

%% for each relevant point on the curve, identify the corresponding point on the boundary
cell_bdry_pt=zeros(numel(s),2);
for si=1:numel(s)
    t1=t(s(si));
    r1t1 = [rtx_test(s(si)) rty_test(s(si))];
    %   Compute r1'(ti) and ||r1'(ti)||
    r1prime_t1 = (4*(1-t1)^3).*(P1-P0)+(12*(1-t1)^2*t1).*(P2-P1)+(12*(1-t1)*t1^2).*(P3-P2)+(4*t1^3).*(P4-P3);
    r1prime_t1_norm=norm(r1prime_t1);
    %   Compute the normal vectors n(ti)
    tangent_t1 = r1prime_t1./(r1prime_t1_norm);
    nt1 = [-tangent_t1(2) tangent_t1(1)];
    
    gamma_val=2*max(max(dist_pair)*10^6);
    R1t1= r1t1+gamma_val.*nt1;
    R1t2= r1t1-gamma_val.*nt1;
    
    l1=[R1t1;R1t2];
    
    [Xi,Yi] = intersections(boundary_pts(:,1),boundary_pts(:,2),l1(:,1),l1(:,2));
    % note: due to the length of the line connecting R1t1 and R1t2, there
    % may be more than one intersection point with the cell boundary. the
    % intersection point that is closest to the point on the curve is the
    % point we want so grab it:
    curve_pts=repmat(r1t1,numel(Xi),1);
    [~,min_idx]=min(sqrt((curve_pts(:,1)-Xi(:,1)).^2+(curve_pts(:,2)-Yi(:,1)).^2));
    cell_bdry_pt(si,:)=[Xi(min_idx(1)), Yi(min_idx(1))];
end

% construct coordinates of all points that represent the boundary of the
% hanging region
Pts_x=[cell_bdry_pt(:,1); rtx_test(s)'];
Pts_y=[cell_bdry_pt(:,2); rty_test(s)'];
Pts_xy=[Pts_x, Pts_y];

% estimate the area of the hanging region (in units um^2)
[~,Ar]=boundary(Pts_xy(:,1),Pts_xy(:,2));

% figure
% plot(Pts_xy(:,1),Pts_xy(:,2),'k*')
% hold on
% plot(Pts_xy(k,1),Pts_xy(k,2),'r')

% outline_scaled=10^6.*outline;
% figure
% line(outline_scaled(1,:),outline_scaled(2,:))
% hold on
% fill(Pts_xy(k,1),Pts_xy(k,2),'r')

% compute energy cost function
E_BD=(1/2)*k_stretch*(Ar/A_FA);
end % end function