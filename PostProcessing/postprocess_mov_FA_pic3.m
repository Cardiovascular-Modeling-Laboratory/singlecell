% File name: postprocess_mov_FA_pic3.m
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
%load(filename,'integrin_bound_store','integrin_free_store','bound_colormap','outline','mat_r','time_store');
load(filename,'integrin_bound_store','integrin_free_store','bound_colormap','outline','mat_r','time_store','outside_segs');

% specify which time point you want to know information for
time_step_u=time_store(2);
thr = input('At what time (in hours) do you want to see the FA distribution? (enter value <=72): ');
tsec=thr*3600; % convert from hours to seconds 
tp=floor(1+tsec/time_step_u); % grab the nearest index which corresponds to the desired time point
time_index = tp;
disp(['Will plot results using the closest timepoint that is smaller than what you asked for: ',num2str(time_store(tp)),' sec'])

filename_pic = [filename '_pic_FA.pdf'];

bound_integrin_temp = integrin_bound_store;
free_integrin_temp = integrin_free_store;
[tl,pl] = size(bound_integrin_temp);

max_bound = max(free_integrin_temp(time_index,:)+bound_integrin_temp(time_index,:));
bound_integrin = 1.*bound_integrin_temp./max_bound;
min_x = min(mat_r(:,1)); %the minimum of the x coordinate
max_x = max(mat_r(:,1)); %the maximum of the x coordinate
min_y = min(mat_r(:,2)); %the minimum of the y coordinate
max_y = max(mat_r(:,2)); %the maximum of the y coordinate
min_z_bound = min(min(bound_integrin));
max_z_bound = max(max(bound_integrin));

clim1 = [0,1];
%clim1 = [min_z_bound,max_z_bound];

h_pl = figure('Visible','off','NextPlot','replacechildren','Color',[1 1 1],'Units','pixels','OuterPosition',[200 200 600 600]);
hold on;
%bound
make_cont_plot_movFA(bound_integrin(time_index,:),mat_r,clim1,{['Normalized bound integrin at t= ',num2str(time_store(tp)),' sec']},outside_segs,outline)
    
%colormap(bound_colormap);
colormap(jet);
colorbar('location','southoutside');
    
tmp = size(outline);
if tmp(1) ==1 %if circle
    rectangle('Position',[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
    %          'FaceColor','r')
elseif tmp(1)==2 % if ellipse
       line(outline(1,:),outline(2,:),'Color',[0 0 0],'LineWidth',2);
    else
    if tmp(2)>5
        rectangle('Position',[outline(1,2),outline(2,2),outline(1,4)-outline(1,3),outline(2,4)-outline(2,2)], 'FaceColor','w','LineStyle','none');
        rectangle('Position',[outline(1,8),outline(2,8),outline(1,6)-outline(1,8),outline(2,6)-outline(2,7)], 'FaceColor','w','LineStyle','none');
    end
    line(outline(1,:),outline(2,:),'Color',[0 0 0],'LineWidth',2);
end

shape = input('Are you running a [1]=square, [2]=stair, [3]=circle, [4]=oval: ');
if shape ==1
   axis([min_x max_x min_y max_y]); %square
%    axis([-0.1 1.0 -0.1 1.0]); %square
elseif shape == 3
    axis([-0.6 0.6 -0.6 0.6]); %circle
elseif shape==4
    axis([-outline(1) outline(1) -outline(1) outline(1)]); 
else
    %axis([-0.01 1.55 -0.1 1.46]); %triangle
%    axis([-0.1 2.0 -0.1 0.85]); %stair
end
hold off;
saveas(h_pl,filename_pic)
close (h_pl)


