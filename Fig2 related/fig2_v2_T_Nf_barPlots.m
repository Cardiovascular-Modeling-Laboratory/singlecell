% this code was used to generate the basic bar plots used in constructing
% figure 2

% load data to be visualized 

meanT_tot=zeros(2,4); %2x4 because there are 2 sets of data to plot on one bar plot, each consistsing of 4 bars
stdT_tot=zeros(2,4);
meanNf_tot=zeros(2,4);
stdNf_tot=zeros(2,4);

T_mat_full=cell(1,2);
Nf_mat_full=cell(1,2);

for i2=1:2
    disp('Please select the results file...');
    [file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\simulation codes');
    filename = [path file];
    load(filename)%,'F_res','dA_store','nsim','AR_vec')%,'F_store','outline','mat_r','dA','Npts_t','time_store','outside_segs','R_store');
    %filename_pic = [filename '_picT_test.pdf'];
    
    T_mat=zeros(size(F_res));
    Nf_mat=zeros(size(F_res));
    for i=1:size(F_res,1)
        for m=1:nsim
            net1=net_res{i,m};
            net=net1{end};
            Nf_vec=zeros(size(net,1),1);
            for j=1:size(net,1)
                net_temp=net{j};
                Nf_vec(j)=size(net_temp,3);
            end
            Nf_mat(i,m)=sum(Nf_vec);
            
            Fs=F_res{i,m};
            T = Fs./dA;
            T_mag = sqrt((T(:,:,1).^2 + T(:,:,2).^2)); % magnitude of traction stress at each point, in N/m^2=Pa
            T_kpa=T_mag.*10^(-3);
            max_z_T = max(max(T_kpa(end,:)));
            T_mat(i,m)=max_z_T;
        end
    end

    % compute the mean and std of the data
    meanT=mean(T_mat,2);
    stdT=zeros(size(meanT));
    meanNf=mean(Nf_mat,2);
    stdNf=zeros(size(meanNf));
    for j=1:numel(stdT)
        stdT(j)=std(T_mat(j,:));
        stdNf(j)=std(Nf_mat(j,:));
    end

    meanT_tot(i2,:)=meanT';
    stdT_tot(i2,:)=stdT';
    meanNf_tot(i2,:)=meanNf';
    stdNf_tot(i2,:)=stdNf'; 
    
    T_mat_full{i2}=T_mat';
    Nf_mat_full{i2}=Nf_mat';

end

%% make bar plot of |T| data
x1=[1 2]; % number of bar plot groups 
mat=colormap(parula(236)); % identify the colors of the bar plots as a matrix of RGB color values

y1=meanT_tot(1,:); % mean data vector 1
X1=repmat(x1',1,numel(y1)); % matrix for grouping data into 2 groups
Y1=nan(size(X1));
Y1(1,:)=y1;
err1=stdT_tot(1,:); % standard deviation vector 1

y2=meanT_tot(2,:); % mean data vector 2
Y1(2,:)=y2;
err2=stdT_tot(2,:); % standard deviation vector 2

% error bar centers need to be aligned relative to the center of the
% bar plots. determine what these centers are
figure
hBar1 = bar(X1,Y1);
ctr_t1=zeros(size(Y1,2),2);
for k1 = 1:size(Y1,2)
    ctr_t1(k1,:) = bsxfun(@plus, hBar1(1).XData, [hBar1(k1).XOffset]');
end

ctr1=ctr_t1(:,1);
ctr2=ctr_t1(:,2);
    
% create plot
figure
bar(X1,Y1);
hold on;
errorbar(ctr1,y1,err1,'k.');
errorbar(ctr2,y2,err2,'k.');
title('|T| for minor, major ')
%%
clear x1 mat y1 y2 Y1 Y2
x1=[1 2]; % number of bar plot groups 
mat=colormap(parula(236)); % identify the colors of the bar plots as a matrix of RGB color values

y1=meanNf_tot(1,:); % mean data vector 1
X1=repmat(x1',1,numel(y1)); % matrix for grouping data into 2 groups
Y1=nan(size(X1));
Y1(1,:)=y1;
err1=stdNf_tot(1,:); % standard deviation vector 1

y2=meanNf_tot(2,:); % mean data vector 2
Y1(2,:)=y2;
err2=stdNf_tot(2,:); % standard deviation vector 2

% error bar centers need to be aligned relative to the center of the
% bar plots. determine what these centers are
figure
hBar1 = bar(X1,Y1);
ctr_t1=zeros(size(Y1,2),2);
for k1 = 1:size(Y1,2)
    ctr_t1(k1,:) = bsxfun(@plus, hBar1(1).XData, [hBar1(k1).XOffset]');
end

ctr1=ctr_t1(:,1);
ctr2=ctr_t1(:,2);
    
% create plot
figure
bar(X1,Y1);
hold on;
errorbar(ctr1,y1,err1,'k.');
errorbar(ctr2,y2,err2,'k.');
title('Nf for minor, major ')

%% statistical analysis for Nf 

U_data=[Nf_mat_full{1} Nf_mat_full{2}]; % data as an array of info where each column corresponds to one group

[p_u,~,stats_u]=anova1(U_data,[],'off'); % the [] and off parts suppress the dispaly of outputs from anova1

if p_u<0.05 % if there's a difference within the group
	% run Tukey-Kramer post-hoc analysis
    bar_compare=multcompare(stats_u,'CType','tukey-kramer','Display','off');
    sig_dif_mat_u=bar_compare(find(bar_compare(:,end)<0.05),[1 2]); % this is an array which shows the groups from U_data that have a significant difference
    p_vals=bar_compare(find(bar_compare(:,end)<0.05),end); % this is an array showing the p-value for each significantly different pair in the array above
    
    % display pairs which are significantly different and the corresponding
    % p-values
    disp('Some pairs have significant differences:')
    disp(sig_dif_mat_u)
    disp(['Corresponding p-vals:  ',num2str(p_vals')])
else
        disp('No significant difference for Nf ')
end

%% statistical analysis for |T| 
clear U_data p_u

U_data=[T_mat_full{1} T_mat_full{2}]; % data as an array of info where each column corresponds to one group

[p_u,~,stats_u]=anova1(U_data,[],'off'); % the [] and off parts suppress the dispaly of outputs from anova1

if p_u<0.05 % if there's a difference within the group
	% run Tukey-Kramer post-hoc analysis
    bar_compare=multcompare(stats_u,'CType','tukey-kramer','Display','off');
    sig_dif_mat_u=bar_compare(find(bar_compare(:,end)<0.05),[1 2]); % this is an array which shows the groups from U_data that have a significant difference
    p_vals=bar_compare(find(bar_compare(:,end)<0.05),end); % this is an array showing the p-value for each significantly different pair in the array above
    
    % display pairs which are significantly different and the corresponding
    % p-values
    disp('Some pairs have significant differences:')
    disp(sig_dif_mat_u)
    disp(['Corresponding p-vals:  ',num2str(p_vals')])
else
        disp('No significant difference for |T|')
end
