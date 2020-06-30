disp('Please select the results file...');
[file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\advanced model+linked');
filename = [path file];
load(filename)%,'F_store','outline','mat_r','dA','Npts_t','time_store','outside_segs','R_store');
%filename_pic = [filename '_picT_test.pdf'];

Lc_cell=cell(size(F_res));
for i=1:size(F_res,1)
    for m=1:nsim
        net1=net_res{i,m};
        net=net1{end};
        Lc_vec=[];
        for j=1:size(net,1)
            net_temp=net{j};
            for k=1:size(net_temp,3)
                net_x_temp=net_temp(1,:,k);
                net_y_temp=net_temp(2,:,k);
                Lc_vec=[Lc_vec, arclength(net_x_temp',net_y_temp')];
            end
        end
        Lc_cell{i,m}=Lc_vec;
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
for j=1:numel(stdLc)
    v1=Lc_cell_combined{j};
    meanLc(j)=mean(v1);
    stdLc(j)=std(v1);
end

x=1:2:11;
figure
errorbar(x,meanLc,stdLc)
xlim([0 12])
xlabel('aspect ratio')
ylabel('Lc in the simulated network')
