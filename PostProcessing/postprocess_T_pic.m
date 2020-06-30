% File name: postprocess_mov.m
% 
% Last updated: 03/29/2010
% First version finished: 03/22/2010
% 
% written by Anya Grosberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Harvard University, Cambridge, MA 02138
% 
% The purpose of this function is to process the results of the single_cell
% simulation and to create output movies and plots
% 
% 
% Input:  1. ask the user to input a file
%        
%
% Global: 1. Fig_Num number of the previous figure number
%        
%         
% Output: countour plot movies


%Declare a variable Fig_Num to be global - this represents the figure
%number, so that each set of figures generates in a different window.
clear
global Fig_Num
Fig_Num = 1;
%
%
disp('Please select the results file...');
[file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\advanced model+linked');
filename = [path file];
load(filename,'F_store','outline','mat_r','dA','Npts_t','time_store','outside_segs','R_store');
filename_pic = [filename '_picT_test.pdf'];

% specify which time point you want to know information for
time_step_u=time_store(2);
thr = input('At what time (in hours) do you want to see the |T| distribution? (enter value <=72): ');
tsec=thr*3600; % convert from hours to seconds 
tp=floor(1+tsec/time_step_u); % grab the nearest index which corresponds to the desired time point
time_index = tp;
disp(['Will plot results using the closest timepoint that is smaller than what you asked for: ',num2str(time_store(tp)),' sec'])


T = F_store./dA; % traction stress in N/m^2=Pa
T_mag = sqrt((T(:,:,1).^2 + T(:,:,2).^2)); % magnitude of traction stress at each point, in N/m^2=Pa
T_kpa=T_mag.*10^(-3); % magnitude of traction stress at each point, in kPa

NT = T_mag/max(max(T_mag)); % this will be used to draw the traction force arrows (later in this code)
min_x = min(mat_r(:,1)); %the minimum of the x coordinate
max_x = max(mat_r(:,1)); %the maximum of the x coordinate
min_y = min(mat_r(:,2)); %the minimum of the y coordinate
max_y = max(mat_r(:,2)); %the maximum of the y coordinate
min_z_T = min(min(T_kpa(time_index,:)));
max_z_T = max(max(T_kpa(time_index,:)));

clim2 = [0,max_z_T];

h_pl = figure('NextPlot','replacechildren','Color',[1 1 1],'Units','pixels','OuterPosition',[200 200 600 600]);
hold on;
%traction force magnitude plot
make_cont_plot_mov(T_kpa(time_index,:),mat_r,clim2,{['Traction Stress (in kPa) at t= ',num2str(time_store(tp)),' sec'],['max|T|=',num2str(max_z_T),' kPa']},outline,outside_segs)
colorbar('location','southoutside');

T_temp_x(:,1) = T(time_index,:,1);
T_temp_y(:,1) = T(time_index,:,2);
%ind = find(NT(length(F_store),:)>1e-2);
ind = find(NT(size(F_store,1),:)>1e-2);
%quiver(mat_r(ind,1),mat_r(ind,2),T_temp_x(ind,1),T_temp_y(ind,1),1,'w','LineWidth',2)
tmp = size(outline);
if tmp(1) ==1 %if circle
    rectangle('Position',[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
    %          'FaceColor','r')
      elseif tmp(1)==2 % if ellipse
       line(outline(1,:),outline(2,:),'Color',[0 0 0],'LineWidth',2);
else
    if tmp(2)>5
        rectangle('Position',[outline(1,2),outline(2,2),outline(1,4)-outline(1,3),outline(2,4)-outline(2,2)], 'FaceColor','w','EdgeColor','w');
        rectangle('Position',[outline(1,8),outline(2,8),outline(1,6)-outline(1,8),outline(2,6)-outline(2,7)], 'FaceColor','w','EdgeColor','w');
    end
    line(outline(1,:),outline(2,:),'Color',[0 0 0],'LineWidth',2);
end
hold off;
saveas(h_pl,filename_pic)
 
    
F_mag=sqrt(F_store(:,:,1).^2+F_store(:,:,2).^2);

r_vals=R_store(tp,:);
lat_sat_idx=intersect(find(F_mag(tp,:)>=3.0e-9), find(r_vals>0.2)); %0.25
%lat_sat_idx=find(F_mag(tp,:)>=3.0e-9);
figure
plot(mat_r(:,1),mat_r(:,2),'k.')
hold on
%plot(mat_r(find(F_mag(tp,:)>=5.5e-9),1),mat_r(find(F_mag(tp,:)>=5.5e-9),2),'bo')
plot(mat_r(lat_sat_idx,1),mat_r(lat_sat_idx,2),'bo')

Nf=zeros(1,size(F_mag,1));
for j=1:size(F_mag,1)
    %Nf(j)=numel(find(F_mag(j,:)>=5.5e-9));
    %Nf(j)=numel(find(F_mag(j,:)>=3.0e-9));
    Nf(j)=numel(intersect(find(F_mag(tp,:)>=3.0e-9), find(R_store(j,:)>0.2)));
end
figure
plot(1:numel(Nf),Nf,'r')
title('Number of saturated points vs time')

    
% % create animation of saturation point evolution
% for k = 1:25:size(F_mag,1)
%     plot(mat_r(:,1),mat_r(:,2),'k.')
%     hold on
%     plot(mat_r(find(F_mag(k,:)>=3.0e-9),1),mat_r(find(F_mag(k,:)>=3.0e-9),2),'bo')
%     title(['t= ',num2str(k),' of ', num2str(size(F_mag,1))])
% 	M(k) = getframe;
%     hold off
% end
% %figure
% %movie(M)