% File name: initial_integrin_units.m
%
% Last updated: 9/11/2018
% Previous versions: 06/24/2010
% First version finished: 08/21/2009
%
% written by Anya Grosberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Harvard University, Cambridge, MA 02138
%
% The purpose of this function is to seed the initial concentration of the
% various integrin types. It does this using the information already
% determined by the simulation. The user will be asked to identify the type
% of initial conditions from a list of possibilities:
%           [1] = Random: All integrins are free and randomly distributed
%           [2] = Uniform: All integrins are free and uniformly distributed
%           [3] = From Previous: to be used if there's already been a
%           simulation performed and user wants the same initial conditions
%           [4] = Random Bound: All integrins are bound to premyofibrils
%           and randomly distributed
%           [5] = Uniform Bound: All integrins are bound to premyofibrils
%           and uniformly distributed
%           [6] = Outline: All integrins are bound to premyofibrils and
%           uniformly distributed along the boundary of the cell
%           [7] = Random Both: Integrins are either free or bound to
%           premyofibrils and are randomly distributed
%           [8] = Uniform Both: Integrins are either free or bound to
%           premyofibrils and are uniformly distributed
%           [9] = Inside Edge Bound: Specific to stair shape cell (ignore)
%           [10] = Random Outline: All integrins are bound to premyofibrils
%           and are  uniformly distributed along the boundary of the cell
%           [11] = Outline Op: Specific to stair shape cell (ignore)
%
%
% Input:  1. Total number of lattice points in the cell, Npts_t (determined
%         by cell_geometry_units file)
%         2. Total amount of integrin, Tot_Integrin_Val (obtained from
%         parameter_gen_units file)
%         3. dA - unit area per point (i.e. A/Npts_t, determined by
%         cell_geometry_units file)
%         4. Concave_ind - a vector the same length as mat_r, but with
%         zeros where there is no ECM in the concave shapes (determined by
%         cell_geometry_units file)
%
%
% Output: 1. Free integrin concentration 
%         2. Pre-myofibril bound integrin concentrations 
%         3. Nascent-myofibril bound integrin cocentration 
%         4. rho_bar_est - estimated average total integrin density
%         determined by this code. This should be approximately equal to
%         the rho_bar value we inputted in parameter_gen_units)
%         5. initial_name - name used to designate the type of IC used

function [integrin_free, integrin_bound, rho_bar, initial_name,percent_err]=initial_integrin_units_v2(Npts_t,Tot_Integrin_Val,dA,Concave_ind,outline,mat_r,rho_bar)

%ask the user of the type of initial conditon they want to run with this
%simulation
initial_choice = input...
    (['What kind of initial condition would you like to start with? 1==[random free integrins], 2==[uniform free integrins], 3==[from previous],'...
    '\n4=[random bound integrins], 5=[uniform bound integrins], 6=[uniform bound integrins on outline], \n7=[random free integinrs and random bound integrins], 8=[uniform free integrins and uniform bound integrins],'...
    '\n9=[bound integrins on inside edge] (FOR STAIR/WEBBED STAIR SHAPE ONLY), 10=[random bound integrins on outline],'... 
    '\n11=[bound integrins on region outside of outline] (FOR WEBBED STAIR SHAPE ONLY): ']);

% CHOICE 3
if initial_choice == 3 %user chose to take the initial conditions from the previous file
    disp('Please select the parameter file that has the initial conditions...');
    [file_ini,path_ini]=uigetfile({'*results_store.mat';'*.*'},'Select Results file with initial conditions File...','G:/MatlabAnalysis/Cell_Model/SingleCell/');
    filename_ini = [path_ini file_ini];
    load(filename_ini,'integrin_free_store','integrin_bound_store');
    integrin_free = integrin_free_store(1,:)';
    integrin_bound = integrin_bound_store(1,:)';
    clear integrin_bound_store integrin_free_store

    initial_name = file_ini(1:(length(file_ini)-18));
% CHOICE 4 OR 5
elseif initial_choice == 4 || initial_choice==5
    %start with a random pre-myofibril integrin configuration
    integrin_free = zeros(Npts_t,1); %no initial pre-myofibril bound integrin
    integrin_n = zeros(Npts_t,1); %no initial nascent myofibril bound integrin

    N_p_cell = Tot_Integrin_Val - dA*sum(integrin_free+integrin_n);

    if initial_choice == 4
        Q = rand(Npts_t,1).*Concave_ind; %Multiply by the concave_ind to get rid of all the points that were outside the ECM island
        initial_name = 'rand_bound';
    end
    if initial_choice==5
        Q= ones(Npts_t,1).*Concave_ind; %Multiply by the concave_ind to get rid of all the points that were outside the ECM island
        initial_name = 'unif_bound';
    end
    sum_Q = sum(Q); %the total now includes only the actual number of points
    integrin_p = (Q./sum_Q).*N_p_cell./dA;
% CHOICE 7 or 8
elseif initial_choice==7 || initial_choice==8
    %have a ratio of bound and unbound integrins
    integrin_n = zeros(Npts_t,1); %no initial nascent myofibril bound integrin

    N_free_cell = rand(1,1).*Tot_Integrin_Val;
    N_p_cell = Tot_Integrin_Val - N_free_cell;

    if initial_choice==7
        Q_free=rand(Npts_t,1).*Concave_ind; %Multiply by the concave_ind to get rid of all the points that were outside the ECM island
        Q_p=rand(Npts_t,1).*Concave_ind;
        initial_name = 'rand_both';
    end
    if initial_choice==8
        Q_free=ones(Npts_t,1).*Concave_ind; %Multiply by the concave_ind to get rid of all the points that were outside the ECM island
        Q_p=ones(Npts_t,1).*Concave_ind;
        initial_name = 'unif_both';
    end
    sum_Q_free = sum(integrin_temp_free);
    sum_Q_p = sum(integrin_temp_p);
    integrin_free = (Q_free./sum_Q_free).*N_free_cell./dA;
    integrin_p = (Q_p./sum_Q_p).*N_p_cell./dA;
% CHOICE 6 or 10
elseif initial_choice==6 || initial_choice==10
    %start with the pre-myofibril integrins at the bounds of the cell
    integrin_free = zeros(Npts_t,1); %no initial pre-myofibril bound integrin
    integrin_n = zeros(Npts_t,1); %no initial nascent myofibril bound integrin
    
    N_p_cell = Tot_Integrin_Val - dA*sum(integrin_free+integrin_n);
    %Set up an array that has ones at the boundary of the stair
    %case and zeros everywhere else
    if length(outline)==5 %square or rectangle cell
        out_ind = ((mat_r(:,1)==outline(1,1))|(mat_r(:,1)==outline(1,2))) | ...
        ((mat_r(:,2)==outline(2,1))|(mat_r(:,2)==outline(2,3)));
    elseif length(outline)==3 % circle cell
        dist_mat=sqrt(mat_r(:,1).^2+mat_r(:,2).^2); % distance between each lattice point ant the center of the circle
        % since we want to isolate the boundary points of the lattice,
        % allow for a little wiggle room (for roundoff error) when
        % determining which distances are "approximately equal" to the
        % circle radius
        out_ind=dist_mat>=0.99.*outline(3);
    elseif length(outline)==4 % isosceles triangle cell
        % here we'll directly compute the lattice indices which are on the
        % boundary
        out_ind=zeros(size(mat_r,1),1);
        yp1=unique(mat_r(:,2)); % isolate the unique y-values
        yp1idx=find(ismember(mat_r(:,2),yp1(1))==1); % grab all indices along the base of the triangle
        out_ind(yp1idx)=1;
        xp2idx=zeros(numel(yp1)-1,2); % for each y-level of the triangle (other than the base), determine the index associated with the min/max x-coord
        for j=2:numel(yp1)
            yp2idx=find(ismember(mat_r(:,2),yp1(j))==1);
            xp2idx(j-1,:)=[min(yp2idx), max(yp2idx)];
        end
        out_ind_temp=unique(xp2idx(:)); % remove duplicate indices
        out_ind(out_ind_temp)=1;
    elseif length(outline)==100 % oval cell
        outline_in=0.99.*outline; % inner oval (oval just slightly inside of the outline)
        outline_out=1.01.*outline; % outter oval (oval just slightly outside of the outline)
        out_ind=zeros(size(mat_r,1),1);
        dist_mat=sqrt(mat_r(:,1).^2+mat_r(:,2).^2); % distance of each lattice point to center of ellipse
        % for each lattice point, construct a ray pointing from the oval's
        % center to a point outside the outter oval. then, determine where
        % this ray interests with the inner and outer ovals. compute the
        % distances of these 2 points. if the lattice point distance is
        % between these 2 computed distances, then it is on the boundary:
        for j=1:size(mat_r,1)
            if sign(mat_r(j,1))==0 % determine if the lattice point is along the vertical axis
                % construct ray
                lx=zeros(1,100);
                ly=linspace(0,3*(0.5+sign(mat_r(j,2)))*max(mat_r(:,2)));
            else
                % determine which quadrant the point is in
                X=mat_r(j,1)>=0;
                Y=mat_r(j,2)>=0;
                quad_number=3+X-Y-2*X*Y;
                % depending on the quadrant, impose the proper sign
                switch quad_number
                    case 1 
                        c=1;
                    case 2 
                        c=-1;
                    case 3
                        c=-1;
                    case 4
                        c=1;
                end
                % construct ray
                xmax=2*c*max(mat_r(:,1));
                lx=linspace(0,xmax);
                ly=linspace(0,(mat_r(j,2)/mat_r(j,1))*xmax);
            end
            [x_in, y_in]=intersections(lx,ly,outline_in(1,:),outline_in(2,:)); % determines where ray interestcs the inner ellipse
            [x_out, y_out]=intersections(lx,ly,outline_out(1,:),outline_out(2,:)); % determines where ray interestcs the outer ellipse
            dist_in=sqrt(x_in.^2+y_in.^2); % distance between inner ntersection point and the center
            dist_out=sqrt(x_out.^2+y_out.^2); % distance between outer ntersection point and the center
            out_ind(j)=((dist_mat(j)>= dist_in)&(dist_mat(j)<=dist_out)); % determine if the lattice point distance is within the distance range
        end        
    else % stair shaped cell
        out_ind = ((mat_r(:,2)==outline(2,1)) | (mat_r(:,1)==outline(1,4)) | (mat_r(:,2)==outline(2,5)) | (mat_r(:,1)==outline(1,1)) ...
        | ((mat_r(:,2)<=outline(2,4)) & (mat_r(:,1)==outline(1,2))) |  ((mat_r(:,1)>outline(1,2)) & (mat_r(:,2)==outline(2,4))) ...
        | ((mat_r(:,2)>outline(2,8)) & (mat_r(:,1)==outline(1,6))) |  ((mat_r(:,1)<=outline(1,6)) & (mat_r(:,2)==outline(2,8))));
    end
    if initial_choice==6
        Q = cast(out_ind,'double');
        initial_name = 'outline';
    end
    if initial_choice==10
        Q = rand(Npts_t,1).*cast(out_ind,'double');
        initial_name = 'outline_rand';
    end
    
    sum_Q= sum(Q); %the total now includes only the actual number of poitns
    integrin_p = (Q./sum_Q).*N_p_cell./dA;
% CHOICE 9
elseif initial_choice==9
    %start with the pre-myofibril integrins at the bounds of the
    %cell on the inside edge
    integrin_free = zeros(Npts_t,1); %no initial pre-myofibril bound integrin
    integrin_n = zeros(Npts_t,1); %no initial nascent myofibril bound integrin

    N_p_cell = Tot_Integrin_Val - dA*sum(integrin_free+integrin_n);
    
    %Set up an array that has ones at the boundary of the stair
    %case and zeros everywhere else
    if length(outline)==5 %square cell
        out_ind = ((mat_r(:,1)==outline(1,1))|(mat_r(:,1)==outline(1,2))) | ...
        ((mat_r(:,2)==outline(2,1))|(mat_r(:,2)==outline(2,3)));
    else
        out_ind = (((mat_r(:,2)<=outline(2,4)) & (mat_r(:,1)==outline(1,2))) | ...
        ((mat_r(:,1)>outline(1,2)) & (mat_r(:,2)==outline(2,4))) ...
        | ((mat_r(:,2)>outline(2,8)) & (mat_r(:,1)==outline(1,6))) |  ((mat_r(:,1)<=outline(1,6)) & (mat_r(:,2)==outline(2,8))));
    end
    Q = cast(out_ind,'double');

    initial_name = 'border_unfav';
    sum_Q = sum(Q); %the total now includes only the actual number of poitns
    integrin_p = (Q./sum_Q).*N_p_cell./dA;
% CHOICE 11
elseif initial_choice==11
    %start with the pre-myofibril integrins at the bounds of the
    %cell on all edges except top and bottom
    integrin_free = zeros(Npts_t,1); %no initial pre-myofibril bound integrin
    integrin_n = zeros(Npts_t,1); %no initial nascent myofibril bound integrin
    
    N_p_cell = Tot_Integrin_Val - dA*sum(integrin_free+integrin_n);
    
    %Set up an array that has ones at the boundary of the stair
    %case and zeros everywhere else
    if length(outline)==5 %square cell
        out_ind = ((mat_r(:,1)==outline(1,1))|(mat_r(:,1)==outline(1,2))) | ...
        ((mat_r(:,2)==outline(2,1))|(mat_r(:,2)==outline(2,3)));
    else
        out_ind = ((mat_r(:,1)==outline(1,1))| (mat_r(:,1)==outline(1,4)) |...
        ((mat_r(:,2)<=outline(2,4)) & (mat_r(:,1)==outline(1,2))) | ...
        ((mat_r(:,1)>outline(1,2)) & (mat_r(:,2)==outline(2,4))) ...
        | ((mat_r(:,2)>outline(2,8)) & (mat_r(:,1)==outline(1,6)))...
        |  ((mat_r(:,1)<=outline(1,6)) & (mat_r(:,2)==outline(2,8))));
    end
    Q = rand(Npts_t,1).*cast(out_ind,'double');

    initial_name = 'outline_unfavor';
    sum_Q = sum(Q); %the total now includes only the actual number of poitns
    integrin_p = (Q./sum_Q).*N_p_cell./dA;
% CHOICE 1 or 2
else %perform common duties for the other two possible initial conditions
    integrin_p = zeros(Npts_t,1); %no initial pre-myofibril bound integrin
    integrin_n = zeros(Npts_t,1); %no initial nascent myofibril bound integrin
    %uniform distribution - because this is a descreet system, I assume that
    %all the integrins are distributed between the points.
    N_free_cell = Tot_Integrin_Val - dA*sum(integrin_p+integrin_n);

    if initial_choice==1 %random distribution
        Q = rand(Npts_t,1).*Concave_ind; %Multiply by the concave_ind to get rid of all the points that were outside the ECM island
        initial_name = 'rand';
    else
        Q = ones(Npts_t,1).*Concave_ind;%uniform distribution of free integrin
        initial_name = 'unif';
    end
    sum_Q = sum(Q); %the total now includes only the actual number of poitns
    integrin_free = (Q./sum_Q).*N_free_cell./dA;
end

integrin_bound=integrin_n+integrin_p;

%average total integrin for a discrete system - the total number of points
%is the points in the ECM island (this should be the same value as the
%input parameter rho_bar and we can use this as a check to make sure we've
%done everything correctly)
Npts_for_binding = Npts_t+(sum(Concave_ind-ones(size(Concave_ind))));
rho_bar_est =(sum(integrin_bound+integrin_free)/Npts_for_binding).*Concave_ind;
percent_err=100*abs(rho_bar_est-rho_bar)/rho_bar;
% if percent_err>1
%     disp('ERROR! Initial integrin values give a rho_bar that is too large!')
% else
%     rho_bar=rho_bar_est;
% end
rho_bar=rho_bar_est;
end