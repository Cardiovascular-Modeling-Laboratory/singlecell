disp('Please select the results file...');
[file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\advanced model+linked');
filename = [path file];
load(filename)%,'F_store','outline','mat_r','dA','Npts_t','time_store','outside_segs','R_store');
%filename_pic = [filename '_picT_test.pdf'];

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
    end
end

% compute the mean and std of the data
meanNf=mean(Nf_mat,2);
stdNf=zeros(size(meanNf));
for j=1:numel(stdNf)
    stdNf(j)=std(Nf_mat(j,:));
end

x=1:2:11;
figure
errorbar(x,meanNf,stdNf)
xlabel('aspect ratio')
ylabel('Number of fibers in the simulated network')





