% this code creates a plot showing how T changes with alpha=0 or alpha=1.
% First, select the case when alpha=1, then select when alpha=0.

meanT_tot=zeros(7,2);
stdT_tot=zeros(7,2);
meanNf_tot=zeros(7,2);
stdNf_tot=zeros(7,2);
meanLc_tot=zeros(7,2);
stdLc_tot=zeros(7,2);


for i2=1:2
    disp('Please select the results file...');
    [file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\advanced model+linked');
    filename = [path file];
    load(filename)%,'F_res','dA_store','nsim','AR_vec')%,'F_store','outline','mat_r','dA','Npts_t','time_store','outside_segs','R_store');
    %filename_pic = [filename '_picT_test.pdf'];
    
    Lc_cell=cell(size(F_res));
    T_mat=zeros(size(F_res));
    Nf_mat=zeros(size(F_res));
    for i=1:size(F_res,1)
        %dA=dA_store(i);
        dA=min(dA_store);
        for m=1:nsim
            net1=net_res{i,m};
            net=net1{end};
            Lc_vec=[];
            Nf_vec=zeros(size(net,1),1);
            for j=1:size(net,1)
                net_temp=net{j};
                for k=1:size(net_temp,3)
                    net_x_temp=net_temp(1,:,k);
                    net_y_temp=net_temp(2,:,k);
                    Lc_vec=[Lc_vec, arclength(net_x_temp',net_y_temp')];
                end
                Nf_vec(j)=size(net_temp,3);
            end
            Nf_mat(i,m)=sum(Nf_vec);
            Lc_cell{i,m}=Lc_vec;
            
            Fs=F_res{i,m};
            T = Fs./dA;
            T_mag = sqrt((T(:,:,1).^2 + T(:,:,2).^2)); % magnitude of traction stress at each point, in N/m^2=Pa
            T_kpa=T_mag.*10^(-3);
            max_z_T = max(max(T_kpa(end,:)));
            T_mat(i,m)=max_z_T;
        end
    end
    
    % combine Lc_cell from an MxN cell to an Mx1 cell where each element is the
    % cumulation of all length terms across all simulations for a given AR
    Lc_cell_combined=cell(size(Lc_cell,1),1);
    for j=1:size(Lc_cell,1)
        Lc_temp=[];
        for k=1:size(Lc_cell,2)
            Lc_temp=[Lc_temp, Lc_cell{j,k}];
        end
        Lc_cell_combined{j}=Lc_temp;
    end

    % compute the mean and std of the data
    meanLc=zeros(1,size(Lc_cell_combined,1));
    stdLc=zeros(1,size(Lc_cell_combined,1));
    
    meanT=mean(T_mat,2);
    stdT=zeros(size(meanT));
    meanNf=mean(Nf_mat,2);
    stdNf=zeros(size(meanNf));
    for j=1:numel(stdT)
        stdT(j)=std(T_mat(j,:));
        stdNf(j)=std(Nf_mat(j,:));
        v1=Lc_cell_combined{j};
        meanLc(j)=mean(v1);
        stdLc(j)=std(v1);
    end

    meanT_tot(:,i2)=meanT;
    stdT_tot(:,i2)=stdT;
    meanNf_tot(:,i2)=meanNf;
    stdNf_tot(:,i2)=stdNf; 
    meanLc_tot(:,i2)=meanLc;
    stdLc_tot(:,i2)=stdLc; 
end

x=AR_vec;%1:2:11;
figure
errorbar(x,meanT_tot(:,1),stdT_tot(:,1))
hold on
errorbar(x,meanT_tot(:,2),stdT_tot(:,2))
%ylim([0,4])
xlim([0,max(x)+1])
legend('a=0','a=1')
axis square
title('max|T| versus cell AR')

x=AR_vec;%1:2:11;
figure
errorbar(x,meanNf_tot(:,1),stdNf_tot(:,1))
hold on
errorbar(x,meanNf_tot(:,2),stdNf_tot(:,2))
%ylim([0,4])
xlim([0,max(x)+1])
legend('a=0','a=1')
axis square
title('Number of fibers versus cell AR')

x=AR_vec;%1:2:11;
figure
errorbar(x,meanLc_tot(:,1),stdLc_tot(:,1))
hold on
errorbar(x,meanLc_tot(:,2),stdLc_tot(:,2))
%ylim([0,4])
xlim([0,max(x)+1])
legend('a=0','a=1')
axis square
title('Length of fibers versus cell AR')

%%
% this code produces a plot of the COOP mean/std as a function of aspect
% ratio.
% Note: only perform this analysis for cell aspect ratios 1,3,5,7,8,11,13.
% For consistency, pick the model results first and then the experimental
% results

%close all
clear all

for i2=1:2
    path = uigetdir('C:\Users\GrosbergLab\Documents\Will\Data to Analyze\','Please pick a base directory...');

    [file,path_p]=uigetfile({'*.mat';'*.*'},['Select COOP results file(s) '],...
        path,'MultiSelect','on');

    if iscell(file) %Determine if the user wants to go through one file or multiple
        NumFiles = length(file);
        for i=1:NumFiles
            temp = file{i};
            FileIndex{i} = temp;
        end
    else
        NumFiles = 1;
        FileIndex{1} = file;
    end

    cell_AR=[1 3 5 7 8 11 13];

    COOP_avg=zeros(1,NumFiles);
    COOP_std=zeros(1,NumFiles);
    AR_vec=zeros(1,NumFiles);

    for count=1:NumFiles

        filename_temp_p = FileIndex{count};
        path_and_filename_p = [path_p filename_temp_p];
        filename_p = path_and_filename_p;
        disp(filename_p)

        %load parameters
        load(filename_p,'S','COOP');

        A=triu(COOP,1);
        At=A.';
        m=tril(true(size(At)),-1);
        COOP_vec=At(m).';

        COOP_avg(count)=mean(COOP_vec);
        COOP_std(count)=std(COOP_vec);

    end

    %
    x = cell_AR';
    y = COOP_avg';
    % Define Start points, fit-function and fit curve
    x0 = [1 1]; 
    %n=1;
    fitfun = fittype( @(a,n,x) (x.^n)./(x.^n+a.^n) );
    [fitted_curve,gof] = fit(x,y,fitfun,'StartPoint',x0);

    % Save the coeffiecient values 
    coeffvals = coeffvalues(fitted_curve);
    a=coeffvals(1);
    n=coeffvals(2);
    Ec50(i2)=a;
    n_val(i2)=n;
    
    xm=linspace(min(x),max(x));
    yf=xm.^n./(xm.^n+a.^n);
    p11 = predint(fitted_curve,xm,0.95,'observation','on');

    % Plot results
    figure
    errorbar(cell_AR,COOP_avg,COOP_std,'ks')
    hold on
    plot(xm,yf,'b')
    plot(xm,p11,'m--')
    xlabel('cell AR')
    ylabel('COOP')
    if i2==1
        title({['model cell COOP values for scale S=',num2str(S)],['Hill function with EC50= ',num2str(a),', n= ',num2str(n)]})
    else
        title({['experiment cell COOP values for scale S=',num2str(S)],['Hill function with EC50= ',num2str(a),', n= ',num2str(n)]})
    end
    ylim([0,1.1])
    hold off
end

% make bar plot of EC50 and n values
figure
bar(Ec50)
title({'EC50 of Hill function', 'left bar: model, right bar: experiment'})

figure
bar(n_val)
title({'Hill coefficient','left bar: model, right bar: experiment'})
