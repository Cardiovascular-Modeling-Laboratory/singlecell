% File name: make_cont_plot.m
% 
% Last updated: 03/29/2009
% First version finished: 09/01/2009
% 
% written by Anya Grosberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Harvard University, Cambridge, MA 02138
% 
% The purpose of this function is to make a countour plot over the cell
% area of any variable passed to this function
% 
% 
% Input:  1. The variable to be plotted as a vector
%         2. The r vector field corresponding to the above variable
%         3. The grid that can be used for plotting
%            cell, i.e. the rectangle that the cell is inscribed in.
%            This grid is needed to calculate the spatial derivatives.
%            The variable (:,:,1) = X grid, (:,:,2) = Y grid
%            rows of the output array X are copies of the vector x;
%            columns of the output array Y are copies of the vector y
%        
%
% Global: 1. Fig_Num number of the previous figure number
%        
%         
% Output: countour plot

function make_cont_plot_pic(z_vec,mat_r,clim,title_st,NC,Mesh)
%Declare a variable Fig_Num to be global - this represents the figure
%number, so that each set of figures generates in a different window.
%global Fig_Num
%
%
Num_points_mesh = Mesh;% the number of points sufficient to make a good mesh for plotting
min_x = min(mat_r(:,1)); %the minimum of the x coordinate
max_x = max(mat_r(:,1)); %the maximum of the x coordinate
dx = (max_x-min_x)/(Num_points_mesh); %the spacing in the x direction of the grid 
min_y = min(mat_r(:,2)); %the minimum of the y coordinate
max_y = max(mat_r(:,2)); %the maximum of the y coordinate
dy = (max_y-min_y)/(Num_points_mesh); %the spacing in the y direction of the grid 
[mesh_cell(:,:,1),mesh_cell(:,:,2)] = meshgrid(min_x:dx:max_x,min_y:dy:max_y); %the mesh of a rectangular spcae surrounding the cell
%grid the data onto the new plotting mesh
Z = griddata(mat_r(:,1),mat_r(:,2),z_vec,mesh_cell(:,:,1),mesh_cell(:,:,2));

%
%Fig_Num = Fig_Num +1;
%figure(Fig_Num)
%Make countour plot
contour(mesh_cell(:,:,1),mesh_cell(:,:,2),Z,NC,'Fill','on','LineStyle','none');%,'LevelListMode','manual','LevelList',LevelList);
%axis equal;
set(gca,'DataAspectRatio',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1],'XTick',[],'YTick',[],'Box','off','CLim',clim);
%title(title_st,'Color',[0 0 0],'FontName','Arial','FontSize',18,'FontWeight','bold');
