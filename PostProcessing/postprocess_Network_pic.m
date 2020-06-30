% this code plots the network at specific time points using the results
% % from "test_0_results_store". You should manually load this file! for
% some reason, it keeps giving an error when trying to load the file...

% ask the user to load the model results file
disp('Please select the results file...');
disp('This may take a while to load because of the file size. BE PATIENT...');
[file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\advanced model+linked');
filename = [path file];
load(filename,'actin_network_t','fiber_id1_t','mat_r','nuc_x','nuc_y','time_store',...
    'outline');
time_step_u=time_store(2);
% specify which time point you want to know information for
thr = input('At what time (in whole hours) do you want to see the network? (enter value <=72): ');
tsec=thr*3600; % convert from hours to seconds 
tp=floor(1+tsec/time_step_u); % grab the nearest index which corresponds to the desired time point
% specify shape so that you plot the cell outline correctly (only
% applicable for circle shape)
tmp = input('Is the cell shape a circle? 1 = yes, 0 = no ');
disp(['Will plot results using the closest timepoint that is smaller than what you asked for: ',num2str(time_store(tp)),' sec'])

%% actin network plot setup
net=actin_network_t{tp};
fid=fiber_id1_t{tp};
Lat=10^6.*mat_r; % convert to micrometer scale since that's the scale the network is saved in

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
    
    % identify the starting and end point of each collection of fibers
    P0=[net_x_temp(1), net_y_temp(1)];
    P4=[net_x_temp(end), net_y_temp(end)];
    starting_pt=find(ismember(Lat,P0,'rows')==1);
    ending_pt=find(ismember(Lat,P4,'rows')==1);
    % use the starting/ending points to record the fiber identity 
    fiber_id_temp=fid{starting_pt, ending_pt};
    fiber_id_vec=[fiber_id_vec; fiber_id_temp'];
end

filename_pic = [filename '_pic_fullNetwork.pdf'];

figure
% plot nucleus (in micrometer scale)
patch(10^6.*nuc_x,10^6.*nuc_y,'blue','FaceAlpha',0.3)
hold on
% plot the network
plot(actin_x',actin_y','k')
% plot cell outline (in micrometer scale)
if tmp(1) ==1 %if circle
    rectangle('Position',10^6.*[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
    %          'FaceColor','r')
else
    line(10^6.*outline(1,:),10^6.*outline(2,:),'Color',[0 0 0],'LineWidth',2);
end
axis square
title({['Network at t = ', num2str(time_store(tp)),' sec']},'FontName','Arial','FontSize',18)


%% plot the networks
% all the fiber information is stored using 2 matrices. matlab can plot
% matrices X and Y but it does so by ploting the columns of X with
% columns of Y. So in the plots of the network, make sure to take the
% transpose
cont=input('Would you like to break the network down by fiber type? 1=yes, 0=no: ');
if cont==1
    % isolate which indices indicate a premyofibril and which indicate a
    % nascent myofibril
    pre_idx_actin=find(fiber_id_vec==1);
    nas_idx_actin=find(fiber_id_vec==2);

    % isolate premyofibril actin network
    actin_pre_x=actin_x(pre_idx_actin,:);
    actin_pre_y=actin_y(pre_idx_actin,:);

    % isolate nascent myofibril actin network
    actin_nas_x=actin_x(nas_idx_actin,:);
    actin_nas_y=actin_y(nas_idx_actin,:);

    choice= input('Which network would you like to see? 1 = Premyofibrils only, 2 = Nascent only, 3 = Full network (Pre+Nascent)');

    switch choice
        case 1
            % -------- plot the actin network only
            filename_pic = [filename '_pic_premyos.pdf'];
            figure
            % plot nucleus (in micrometer scale)
            patch(10^6.*nuc_x,10^6.*nuc_y,'blue','FaceAlpha',0.3)
            hold on
            % plot cell outline (in micrometer scale)
            if tmp(1) ==1 %if circle
                rectangle('Position',10^6.*[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
                %          'FaceColor','r')
            else
                line(10^6.*outline(1,:),10^6.*outline(2,:),'Color',[0 0 0],'LineWidth',2);
            end
            axis square
            % plot the premyofibrils in magenta
            plot(actin_pre_x',actin_pre_y','m')
            title({['Premyofibril network at t = ', num2str(time_store(tp)),' sec']},'FontName','Arial','FontSize',18)
        case 2
            % -------- plot the nascent network only
            filename_pic = [filename '_pic_nascentmyos.pdf'];
            figure
            % plot nucleus (in micrometer scale)
            patch(10^6.*nuc_x,10^6.*nuc_y,'blue','FaceAlpha',0.3)
            hold on
            % plot cell outline (in micrometer scale)
            if tmp(1) ==1 %if circle
                rectangle('Position',10^6.*[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
                %          'FaceColor','r')
            else
                line(10^6.*outline(1,:),10^6.*outline(2,:),'Color',[0 0 0],'LineWidth',2);
            end
            axis square
            % plot the nascent myofibrils in green
            plot(actin_nas_x',actin_nas_y','g')
            title({['Nascent myofibril network at t = ', num2str(time_store(tp)),' sec']},'FontName','Arial','FontSize',18)
        case 3
            % -------- plot the full actin network
            filename_pic = [filename '_pic_fullNetwork.pdf'];

            figure
            % plot nucleus (in micrometer scale)
            patch(10^6.*nuc_x,10^6.*nuc_y,'blue','FaceAlpha',0.3)
            hold on
            % plot cell outline (in micrometer scale)
            if tmp(1) ==1 %if circle
                rectangle('Position',10^6.*[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
                %          'FaceColor','r')
            else
                line(10^6.*outline(1,:),10^6.*outline(2,:),'Color',[0 0 0],'LineWidth',2);
            end
            axis square
            % plot the actin network: premyofibrils in magenta, nascent myofibrils in
            % green
            if ~isempty(actin_pre_x) && isempty(actin_nas_x)
                plot(actin_pre_x',actin_pre_y','m')
            elseif ~isempty(actin_nas_x) && isempty(actin_pre_x)
                plot(actin_nas_x',actin_nas_y','g')
            elseif ~isempty(actin_pre_x) && ~isempty(actin_nas_x)
                plot(actin_pre_x',actin_pre_y','m',actin_nas_x',actin_nas_y','g')
            end
            axis([min(Lat(:,1)) max(Lat(:,1)) min(Lat(:,1)) max(Lat(:,1))])
            axis square
            title({['Full network network at t = ', num2str(time_store(tp)),' sec'],['\color{magenta}Premyofibril \color{green} Nascent']},'FontName','Arial','FontSize',18)
            %        title(['Full network network at t = ', num2str(time_store(tp)),' sec'],'FontName','Arial','FontSize',18)
        otherwise
            disp('You picked a bad value. Try again!')
    end
end


