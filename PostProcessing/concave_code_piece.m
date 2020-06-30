[mesh_cell(:,:,1),mesh_cell(:,:,2)] = meshgrid(min_x:dx:max_x,min_y:dy:max_y); %the mesh of a rectangular spcae surrounding the cell
X=mesh_cell(:,:,1);
Y=mesh_cell(:,:,2);
%% test code
x1=unique(mat_r(:,1)); y1=unique(mat_r(:,2));
[X1,Y1] = meshgrid(x1,y1);
k=inpoly([X1(:),Y1(:)],outline');
k2=find(k==0);
X2=[mat_r(:,1); X1(k2)];
Y2=[mat_r(:,2); Y1(k2)];
Z2=[z_vec'; nan(size(X1(k2)))];
Z=griddata(X2,Y2,Z2,X,Y);

figure
%Make countour plot
contour(mesh_cell(:,:,1),mesh_cell(:,:,2),Z,'Fill','on','LineStyle','none');%,'LevelListMode','manual','LevelList',LevelList);
hold on
line(outline(1,:),outline(2,:))
plot(X2,Y2,'ro',X1(k2),Y1(k2),'b*')

%%
lx=min_x:dx:max_x;
ly=min_y:dy:max_y;
% include boundary points so that they'll show up on the meshed grid
lx=lx;
ly=ly;
[mesh_cell(:,:,1),mesh_cell(:,:,2)] = meshgrid(lx,ly); %the mesh of a rectangular spcae surrounding the cell
X=mesh_cell(:,:,1);
Y=mesh_cell(:,:,2);
Z=griddata(mat_r(:,1),mat_r(:,2),z_vec,X,Y);

% % for concave shapes, make the z-values NaN for all points that lie outside
% % the cell boundary:
k=inpoly([X(:),Y(:)],outline');
k2=find(k==0);
Z(k2)=nan;

figure
%Make countour plot
contour(mesh_cell(:,:,1),mesh_cell(:,:,2),Z,'Fill','on','LineStyle','none');%,'LevelListMode','manual','LevelList',LevelList);
set(gca,'DataAspectRatio',[1 1 1],'XColor',[1 1 1],'YColor',[1 1 1],'XTick',[],'YTick',[],'Box','off','CLim',[0 max(z_vec)]);
hold on
line(outline(1,:),outline(2,:))
%plot(X2,Y2,'ro',X1(k2),Y1(k2),'b*')
plot(X,Y,'k.')
%%
