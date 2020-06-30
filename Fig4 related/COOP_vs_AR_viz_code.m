% this code creates the coop plot as a function of aspect ratio. For each
% case (model-model, model-exp, exp-exp), identify the small (S=6) and
% large (S=93) coop results files. This code will then creeate a plot of
% the coop for each case with small scale in red and large scale in blue

clear all

for i2=1:3
    for j2=1:2
        if i2==1
            path = uigetdir('C:\Users\GrosbergLab\Documents\Will\Data to Analyze\','Please pick a base directory for MODEL-MODEL comparisons...');
        elseif i2==2
            path = uigetdir('C:\Users\GrosbergLab\Documents\Will\Data to Analyze\','Please pick a base directory for EXPERIMENT-EXPERIMENT comparisons...');
        else
            path = uigetdir('C:\Users\GrosbergLab\Documents\Will\Data to Analyze\','Please pick a base directory for MODEL-EXPERIMENT comparisons...');
        end

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

        if j2==1
            COOP_avg=zeros(2,NumFiles);
            COOP_std=zeros(2,NumFiles);
            AR_vec=zeros(1,NumFiles);
        end

        for count=1:NumFiles
            Ar=cell_AR(count);

            filename_temp_p = FileIndex{count};
            path_and_filename_p = [path_p filename_temp_p];
            filename_p = path_and_filename_p;
            disp(filename_p)

            %load parameters
            load(filename_p,'S','COOP');

            if i2==3
                % the COOP matrix contains all possible combinations (model-model,
                % model-exp, and exp-exp). isolate the model-exp submatrices
                i_me=1:5;
                j_me=6:10;

                COOP_isolated_me=COOP(i_me,j_me); % COOP matrix containing only COOP values for model-experiment pairs
                C_me=COOP_isolated_me;

                % record COOP avg/std for model-exp comparisions
                COOP_vec=C_me(:);
            else
                A=triu(COOP,1);
                At=A.';
                m=tril(true(size(At)),-1);
                COOP_vec=At(m).';
            end


            COOP_avg(j2,count)=mean(COOP_vec);
            COOP_std(j2,count)=std(COOP_vec);
        end
    
        % Now determine fitted curve with 95% confidence interval
        x = cell_AR';
        y = COOP_avg(j2,:)';
        % Define Start points, fit-function and fit curve
        
        if j==1
            x0 = [1 1]; 
        else
            x0=[1 5];
        end
        
        %n=1;
        fitfun = fittype( @(a,n,x) (x.^n)./(x.^n+a.^n) );
        [fitted_curve,gof] = fit(x,y,fitfun,'StartPoint',x0);

        % Save the coeffiecient values 
        coeffvals = coeffvalues(fitted_curve);
        a=coeffvals(1);
        n=coeffvals(2);
    
     xm=linspace(min(x),max(x));
     yf(j2,:)=xm.^n./(xm.^n+a.^n);
%     p11 = predint(fitted_curve,xm,0.95,'observation','on');
    end
    
    % Plot results
    figure
    for j2=1:2
        errorbar(cell_AR,COOP_avg(j2,:),COOP_std(j2,:),'k.') % plot COOP data
        hold on
        plot(xm,yf) % plot fitted curve
    %    plot(xm,p11,'m--') % plot confidence bounds
        xlabel('cell AR')
        ylabel('COOP')
    end
    if i2==1
        title({['model-model cell COOP values for scale S=',num2str(S)],['Hill function with EC50= ',num2str(a),', n= ',num2str(n)]})
    elseif i2==2
        title({['experiment-experiment COOP values for scale S=',num2str(S)],['Hill function with EC50= ',num2str(a),', n= ',num2str(n)]})
    else
        title({['model-experiment COOP values for scale S=',num2str(S)],['Hill function with EC50= ',num2str(a),', n= ',num2str(n)]})
    end
    ylim([0,1.1])
    hold off


end

