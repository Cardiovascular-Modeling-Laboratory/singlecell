% for each aspect ratio, the full coop matrix contains all pairs (exp-exp,
% model-model, model-exp, exp-model). Since model-exp is the same as
% exp-model and that's what we're interested in, we only need a sub-block
% of the full matrix. For each aspect ratio, we've identified the needed
% subblock where each row is the model sim and each column is the exp cell
% Ar = COOP is 10x10, need block (1:5,6:10)

cell_AR=[1 3 5 7 8 11 13];
[file,path_p]=uigetfile({'*.mat';'*.*'},['Select COOP .mat result files for MODEL-EXPERIMENT  '],...
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

% create storage for the average/std COOP values of the model/exp
% comparisions. Each row will correspdond to an aspect ratio. The columns
% indicate which data pairing is being considered: column 1=model-model,
% column 2=model-exp, column 3=exp-exp
COOP_avg_pairs=zeros(NumFiles,3);
COOP_std_pairs=zeros(NumFiles,3);

for count=1:NumFiles
    Ar=cell_AR(count);
    
    filename_temp_p = FileIndex{count};
    path_and_filename_p = [path_p filename_temp_p];
    filename_p = path_and_filename_p;
    
    %load parameters
    load(filename_p,'S','COOP')

    % Produce the input matrix data
    C_full=COOP;

    % Set [min,max] value of C to scale colors
    clrLim = [0,1]; 

    % Set the  [min,max] of diameter where 1 consumes entire grid square
    diamLim = [1, 1];

    % Compute center of each circle 
    % This assumes the x and y values were not entered in imagesc()
    x = 1 : 1 : size(C_full,2); % x edges
    y = 1 : 1 : size(C_full,1); % y edges
    [xAll, yAll] = meshgrid(x,y); 

    % Set color of each rectangle
    cmap = parula(256); 
    Cscaled = (C_full - clrLim(1))/range(clrLim); % always [0:1]
    colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));

    % Set size of each circle
    diamSize = Cscaled * range(diamLim) + diamLim(1); 
%     diamSize = triu(Cscaled * range(diamLim) + diamLim(1)); % uncomment this instead of above line if using C_mm or C_ee

    % Create figure
    fh = figure(); 
    ax = axes(fh); 
    hold(ax,'on')
    colormap(ax,'parula');

    % Create circles
    theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
    h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
        diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:)),1:numel(xAll)); 
    axis(ax,'equal')
    axis(ax,'tight')
    set(ax,'YDir','Reverse')
    colorbar()
    caxis(clrLim);
    title({'COOP values for Model (row) vs Experiment (column)',['AR=',num2str(Ar),', S=',num2str(S)]})
    
    % next, create the bar plot figure showing the mean COOP value for
    % model-model comparisions, model-exp comparisions, and exp-exp
    % comparisions

    i_me=1:5;
    j_me=6:10;
    i_mm=1:5;
    j_mm=1:5;
    i_ee=6:10;
    j_ee=6:10;
    
    COOP_isolated_me=COOP(i_me,j_me); % COOP matrix containing only COOP values for model-experiment pairs
    COOP_isolated_mm=COOP(i_mm,j_mm); % COOP matrix containing only COOP values for model-model pairs
    COOP_isolated_ee=COOP(i_ee,j_ee); % COOP matrix containing only COOP values for experiment-experiment pairs
    
    C_me=COOP_isolated_me;
    C_mm=triu(COOP_isolated_mm,1); % isolate upper triangular portion because this matrix is symmetric with 1s along diagonal
    C_ee=triu(COOP_isolated_ee,1); % isolate upper triangular portion because this matrix is symmetric with 1s along diagonal
    
    % record COOP avg/std for model-model comparisions
    At=C_mm.';
    m=tril(true(size(At)),-1);
    COOP_vec_mm=At(m).';
    COOP_avg_pairs(count,1)=mean(COOP_vec_mm);
    COOP_std_pairs(count,1)=std(COOP_vec_mm);

    % record COOP avg/std for model-exp comparisions
    COOP_avg_pairs(count,2)=mean(C_me(:));
    COOP_std_pairs(count,2)=std(C_me(:));
    
    % record COOP avg/std for exp-exp comparisions
    At=C_ee.';
    m=tril(true(size(At)),-1);
    COOP_vec_ee=At(m).';
    COOP_avg_pairs(count,3)=mean(COOP_vec_ee);
    COOP_std_pairs(count,3)=std(COOP_vec_ee);
    
    figure
    bar(COOP_avg_pairs(count,:))
    barnames={'model-model'; 'model-exp'; 'exp-exp' };
    set(gca,'xticklabel',barnames)
    title({'COOP values for Model/Experiment comparisons',['AR=',num2str(Ar),', S=',num2str(S)]})
    ylim([0,1])
    
    % run anova statistical analysis using grouped data
    data_set_full = [C_me(:)' COOP_vec_ee];
    data_grouping=[ones(size(C_me(:)')) 2.*ones(size(COOP_vec_ee))  ];
    
    [p_u,~,stats_u]=anova1(data_set_full,data_grouping,'off');
    
    if p_u<0.05 % if there's a difference somewhere
        % run Tukey-Kramer post-hoc analysis
        bar_compare=multcompare(stats_u,'CType','tukey-kramer','Display','off');
        sig_dif_mat_u=bar_compare(find(bar_compare(:,end)<0.05),[1 2]); % this is an array which shows the groups from U_data that have a significant difference
        p_vals=bar_compare(find(bar_compare(:,end)<0.05),end); % this is an array showing the p-value for each significantly different pair in the array above

        % display pairs which are significantly different and the corresponding
        % p-values
        disp(['Some pairs have significant differences for AR=',num2str(Ar)])
        disp(sig_dif_mat_u)
        disp(['Corresponding p-vals:  ',num2str(p_vals')])
    else
            disp(['No significant difference for AR=',num2str(Ar)])
    end

end

%% make a colorbar plot of mean coop values

    % Produce the input matrix data
    C_full=COOP_avg_pairs;

    % Set [min,max] value of C to scale colors
    clrLim = [0,1]; 

    % Set the  [min,max] of diameter where 1 consumes entire grid square
    diamLim = [min(min(COOP_std_pairs)), max(max(COOP_std_pairs))];

    % Compute center of each circle 
    % This assumes the x and y values were not entered in imagesc()
    x = 1 : 1 : size(C_full,2); % x edges
    y = 1 : 1 : size(C_full,1); % y edges
    [xAll, yAll] = meshgrid(x,y); 

    % Set color of each rectangle
    cmap = parula(256); 
    Cscaled = (C_full - clrLim(1))/range(clrLim); % always [0:1]
    colIdx = discretize(Cscaled,linspace(0,1,size(cmap,1)));

    % Set size of each circle
    diamSize = Cscaled * range(diamLim) + diamLim(1); 
 % Create figure
    fh = figure(); 
    ax = axes(fh); 
    hold(ax,'on')
    colormap(ax,'parula');

    % Create circles
    theta = linspace(0,2*pi,50); % the smaller, the less memory req'd.
    h = arrayfun(@(i)fill(diamSize(i)/2 * cos(theta) + xAll(i), ...
        diamSize(i)/2 * sin(theta) + yAll(i), cmap(colIdx(i),:)),1:numel(xAll)); 
    axis(ax,'equal')
    axis(ax,'tight')
    set(ax,'YDir','Reverse')
    colorbar()
    caxis(clrLim);
    barnames={'model-model'; 'model-exp'; 'exp-exp' };
    set(gca,'xticklabel',barnames)
%    title({'COOP values for Model-model (row) vs Experiment (column)',['AR=',num2str(Ar),', S=',num2str(S)]})

