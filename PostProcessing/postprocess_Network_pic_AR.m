% this code creates the fiber network plot using the results from the
% AR simulation. 

t=linspace(0,1);
j2=2;
scale_x=linspace(0,1).*10^(-5);
scale_y=ones(size(scale_x)).*(4.9e-5);

outline=outline_store{j2,1};
for m=1:nsim %m3=1:numel(n1)
%    m=n1(m3);%m2(m3);
    nuc_cx=nuc_cx_store(j2,m);
    nuc_cy=nuc_cy_store(j2,m);
    % nucleus curve
    nuc_x=nuc_radius.*cos(2.*pi.*t)+nuc_cx;
    nuc_y=nuc_radius.*sin(2.*pi.*t)+nuc_cy;

    net1=net_res{j2,m};
    net=net1{end};
    % convert actin network into 2 matrix arrays, containing x coords and y
    % coords
    actin_x=[]; % matrix that will hold xcoords of every fiber
    actin_y=[]; % matrix that will hold ycoords of every fiber
    fiber_id_vec=[]; % matrix that will hold fiber identification of every fiber
    for j=1:size(net,1)
        net_temp=net{j};
        for k=1:size(net_temp,3)
            net_x_temp=net_temp(1,:,k);
            net_y_temp=net_temp(2,:,k);
            actin_x=[actin_x; net_x_temp];
            actin_y=[actin_y; net_y_temp];
        end
    end
    
    figure
    tmp = size(outline);
    hold on;
    if tmp(1) ==1 %if circle; circle behaves weird with "patch" so make nucleus a filled circle
        rectangle('Position',[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
        rectangle('Position',[nuc_cx-nuc_radius,nuc_cy-nuc_radius,2*nuc_radius,2*nuc_radius],'Curvature',[1,1],'FaceColor','b');% this line plots the nucleus as a filled blue circle
    else
         line(outline(1,:),outline(2,:),'Color',[0 0 0],'LineWidth',2);
%        patch(nuc_x,nuc_y,'blue','facealpha',0.3); % this line plots the nucleus as a transparant blue circle
    end
    hold on
        plot(10^(-6).*actin_x',10^(-6).*actin_y','k','LineWidth',1)
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        plot(scale_x,scale_y,'r','LineWidth',2) % 10um scale bar
%     xlim([0 50].*10^(-6))
%     ylim([0 50].*10^(-6))
        hold off;
    axis equal;
    set(gca,'Units','centimeters','OuterPosition',[5,5,6,6],'Position',[3,1,9.5,9.5]);
    hold off
    %axis square
end

