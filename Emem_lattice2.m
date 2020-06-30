% this code determines which lattice points the proposed curve passes near
% and the resulting membrane deformation energy 
function [E_mem,lattice_sat_vals_temp]=Emem_lattice2(rtx,rty,Lat,dA_obstruct_radius,df,lattice_sat_vals)

Npts=size(Lat,1); % number of lattice points in the discretized cell
r_col_x = repmat(Lat(:,1),1,Npts); %make the array with each column having rx coorinates
r_row_x = r_col_x'; %make the array with each row having rx coorinates
r_col_y = repmat(Lat(:,2),1,Npts); %make the array with each column having ry coorinates
r_row_y = r_col_y';
drx = r_row_x - r_col_x; %the differences in x coordinates - each row - vary r', column - vary r
dry = r_row_y - r_col_y;
diag_indices=linspace(1,Npts^2,Npts); 
drx(diag_indices)=2*max(max(drx));
dry(diag_indices)=2*max(max(dry));

% necessary parameters
k_b=20; % ~20k_B*T; on order of tens
sig=3.5/10^(-6); %~ 3.5k_B*T/nm^2 converted to k_B*T/um^2; can range from 0 up to 3.5 k_B*T/nm^2
dx_step=min(min(drx));
dy_step=min(min(dry));

% determine which lattice points the curve is near, i.e., those lattice
% points for which the curve is within dA_obstruct_radius 
curve(:,1)=rtx;
curve(:,2)=rty;

[~,distances] = dsearchn(curve,Lat); % this built in function determines the distances between the points in Lat and the proposed curve

idx=find(distances<=dA_obstruct_radius); % find lattice indices that are close to the curve

if isempty(idx)
    figure
    plot(Lat(:,1),Lat(:,2),'r*')
    figure
    plot(curve(:,1),curve(:,2),'b')
end

% increment the fiber number value of those lattice points by 1
% (temporarily because we don't know if we're going to keep this curve just
% yet)
lattice_sat_vals_temp=lattice_sat_vals;
lattice_sat_vals_temp(idx)=lattice_sat_vals_temp(idx)+1;

% record height data for the current membrane energy and the proposed
% membrane energy
h_temp_vec=df.*lattice_sat_vals_temp;
h_vec=df.*lattice_sat_vals;

% in order to compute the spatial derivatives, we need to turn the h
% functions into square matrices
[X,Y]=meshgrid(unique(Lat(:,1)),unique(Lat(:,2)));
h_temp=griddata(Lat(:,1),Lat(:,2),h_temp_vec,X,Y);
h=griddata(Lat(:,1),Lat(:,2),h_vec,X,Y);

% in the convex shapes, the meshgrid will create points that don't have any
% information. these values may be set to NaN in griddata. Set these values
% to 0 
h(isnan(h))=0;
h_temp(isnan(h_temp))=0;

% compute the gradient of h (h_x,h_y) for each h function
[h_x_temp, h_y_temp] = gradient(h_temp,dx_step,dy_step);
[h_x, h_y] = gradient(h,dx_step,dy_step);

% compute the second derivative (h_xx, h_yy) for each h function
h_xx_temp=gradient(h_x_temp,dx_step);
h_yy_temp=gradient(h_y_temp,dy_step);
h_xx=gradient(h_x,dx_step);
h_yy=gradient(h_y,dy_step);

% compute E_mem at each x for both h functions
h_xx_temp_vec=zeros(1,Npts);
h_yy_temp_vec=zeros(1,Npts);
h_x_temp_vec=zeros(1,Npts);
h_y_temp_vec=zeros(1,Npts);

h_xx_vec=zeros(1,Npts);
h_yy_vec=zeros(1,Npts);
h_x_vec=zeros(1,Npts);
h_y_vec=zeros(1,Npts);

for j=1:Npts
    x0 = Lat(j,1);
    y0 = Lat(j,2);
    t = (X == x0) & (Y == y0);
    indt = find(t);

    h_xx_temp_vec(j) = h_xx_temp(indt);
    h_yy_temp_vec(j) = h_yy_temp(indt);
    h_x_temp_vec(j) = h_x_temp(indt);
    h_y_temp_vec(j) = h_y_temp(indt);
    
    h_xx_vec(j) = h_xx(indt);
    h_yy_vec(j) = h_yy(indt);
    h_x_vec(j) = h_x(indt);
    h_y_vec(j) = h_y(indt);
end

E_mem_x_temp=(k_b/2).*(h_xx_temp_vec+h_yy_temp_vec).^2+(sig/2).*(h_x_temp_vec.^2+h_y_temp_vec.^2);
E_mem_x=(k_b/2).*(h_xx_vec+h_yy_vec).^2+(sig/2).*(h_x_vec.^2+h_y_vec.^2);

% compute E_mem
dE_mem_x=abs(E_mem_x_temp-E_mem_x);
E_mem=sum(dE_mem_x)./Npts; % E_mem =(1/A)*sum(dE_mem_x)*dA but dA=A/Npts so this simplifies to the expression in the code

end % end function handle


