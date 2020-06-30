disp('Please select the results file...');
[file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\advanced model+linked');
filename = [path file];
load(filename,'F_res','dA_store','nsim','AR_vec')%,'F_store','outline','mat_r','dA','Npts_t','time_store','outside_segs','R_store');
%filename_pic = [filename '_picT_test.pdf'];

T_mat=zeros(size(F_res));
for i=1:size(F_res,1)
    %dA=dA_store(i);
    dA=min(dA_store);
    for m=1:nsim
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
for j=1:numel(stdT)
    stdT(j)=std(T_mat(j,:));
end

x=AR_vec;%1:2:11;
figure
errorbar(x,meanT,stdT)
ylim([0,4])
xlim([0,max(x)+1])
axis square





