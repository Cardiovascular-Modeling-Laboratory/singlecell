% this code is designed to show how the R functions look at the time point
% of interest.
disp('Please select the results file...');
[file,path]=uigetfile({'*results_store.mat';'*.*'},'Select cell results file...','C:\Users\Grosberg CML\Documents\Will\advanced model+linked');
filename = [path file];
load(filename,'R_store','Npts_t','time_store','Rsat_lat');
filename_pic = [filename '_picR_test.pdf'];

% specify which time point you want to know information for
time_step_u=time_store(2);
thr = input('At what time (in hours) do you want to see the FA distribution? (enter value <=72): ');
tsec=thr*3600; % convert from hours to seconds 
tp=floor(1+tsec/time_step_u); % grab the nearest index which corresponds to the desired time point
time_index = tp;
disp(['Will plot results using the closest timepoint that is smaller than what you asked for: ',num2str(time_store(tp)),' sec'])

Rs=R_store(time_index,:);
pts=1:Npts_t;

rsatq=input('Is Rsat_lat in this file? yes=1, no=0 ');

if rsatq==1
    figure
    plot(pts,Rs,'k-o',pts,Rsat_lat.*ones(size(pts)),'g--')
    xlabel('lattice points')
    ylabel('R function values')
    legend('R(x)','Rsat_{lat}')
    title(['R fcn plots at ',num2str(thr),' hrs'])
else
    figure
    plot(pts,Rs,'k-o')
    xlabel('lattice points')
    ylabel('R function values')
    legend('R(x)')
    title(['R fcn plots at ',num2str(thr),' hrs'])
end
