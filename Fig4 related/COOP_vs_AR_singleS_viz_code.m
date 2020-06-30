% this code creates the coop plot as a function of aspect ratio for each
% case (model-model, model-exp, exp-exp) at a single coop scale. for the
% paper, I used the S=18 scale. note: this isn't plotting the actual data
% but a Hill function fit of the data

clear all

for sval=1:3
    COOP_fit=zeros(3,100);
    for i2=1:3
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


            COOP_avg(1,count)=mean(COOP_vec);
        end
    
        % Now determine fitted curve with 95% confidence interval
        x = cell_AR';
        y = COOP_avg(1,:)';
        % Define Start points, fit-function and fit curve
        if S==93
            x0=[1 5];
        else
            x0 = [1 1]; 
        end
               
        fitfun = fittype( @(a,n,x) (x.^n)./(x.^n+a.^n) );
        [fitted_curve,gof] = fit(x,y,fitfun,'StartPoint',x0);

        % Save the coeffiecient values 
        coeffvals = coeffvalues(fitted_curve);
        a=coeffvals(1);
        n=coeffvals(2);
        
        xm=linspace(min(x),max(x));
        COOP_fit(i2,:)=xm.^n./(xm.^n+a.^n);
        if i2==2
            p11 = predint(fitted_curve,xm,0.95);%,'observation','on');
        end
    end
    
    % Plot results
    figure
    x2 = [xm, fliplr(xm)];
    inBetween = [p11(:,1)', fliplr(p11(:,2)')];
    fill(x2, inBetween, 'g');
    hold on
    for j2=1:3
        plot(xm,COOP_fit(j2,:)) % plot fitted curve
        xlabel('cell AR')
        ylabel('COOP')
    end
    ylim([0,1.1])
    legend({'MM','EE','ME'})
%    plot(xm,p11,'m--') % plot confidence bounds
    title(['COOP values for scale S=',num2str(S)])
    hold off

end

 figure
    hold on
    for j2=1:3
        plot(xm,COOP_fit(j2,:)) % plot fitted curve
        xlabel('cell AR')
        ylabel('COOP')
    end
    ylim([0,1.1])
    legend({'MM','EE','ME'})
%    plot(xm,p11,'m--') % plot confidence bounds
    title(['COOP values for scale S=',num2str(S)])
    hold off
