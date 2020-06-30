disp('Please select the results file...');
[file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\advanced model+linked');
filename = [path file];
load(filename,'combo_res','dA_store','nsim','A_vec')%,'F_store','outline','mat_r','dA','Npts_t','time_store','outside_segs','R_store');
%filename_pic = [filename '_picT_test.pdf'];

Nv=zeros(size(combo_res));
for k=1:6
    for j=1:6
       Nv(k,j)= size(combo_res{k,j}{1,end},1).*dA_store(k).*10^12;
    end
end

Nv_mean=zeros(size(Nv,1),1);
Nv_std=zeros(size(Nv,1),1);

for k=1:6
    Nv_mean(k)=mean(Nv(k,:));
    Nv_std(k)=std(Nv(k,:));
end

figure
errorbar(A_vec.*10^12,Nv_mean,Nv_std)
xlabel('square cell area (um)')
ylabel('estimated focal adhesion area (um)')