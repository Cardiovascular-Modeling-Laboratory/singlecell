% this code runs single_cell_units_linked_v3 multiple times on
% square/rectangle cells with the nucleus placement randomized.
% in the resuls, each row corresonds to one aspect ratio (1,3,5,7,9, or 11)
% and each column is one simulation for a given aspect ratio

clear;

nsim=6; % max number of simulations
choice_vec=[1 2 2 2 2 2 2]; % designate the cell shapes
AR_vec=[1 3 5 7 8 11 13]; % designate the aspect ratios
nuc_decide_vec=[2 2 2 2 2 2 2]; % vector indicating if the nucleus should be placed in the center (1), randomly (2), or manually (3)
nuc_cx_temp_vec=[0 0 0 0 0 0 0]; % vector indicating how x-coord of nucleus should be offset if being placed manually
nuc_cy_temp_vec=[0 0 0 0 0 0 0]; % vector indicating how y-coord of nucleus should be offset if being placed manually
nuc_rel_vec=[1 1 1 1 1 1 1]; % vector indicating if nucleus is treated as an obstruction (1) or not (0)

dA_store=zeros(size(AR_vec));

net_res=cell(numel(nuc_rel_vec),nsim);
F_res=cell(numel(nuc_rel_vec),nsim); 
fiber_id_res=cell(numel(nuc_rel_vec),nsim); 
combo_res=cell(numel(nuc_rel_vec),nsim);
outline_store=cell(numel(nuc_rel_vec),nsim);
nuc_cx_store=zeros(numel(nuc_rel_vec),nsim);
nuc_cy_store=zeros(numel(nuc_rel_vec),nsim);


tic

for my2=1:numel(nuc_rel_vec)
    disp(['AR number ',num2str(my2)])
    nuc_rel=nuc_rel_vec(my2);
    nuc_decide=nuc_decide_vec(my2);
    nuc_cx_temp=nuc_cx_temp_vec(my2);
    nuc_cy_temp=nuc_cy_temp_vec(my2);
    choice=choice_vec(my2);
    ARtemp=AR_vec(my2);
    
    for my1=1:nsim
        disp(['  sim number ',num2str(my1)])
        %load parameters
        load('paramTEST5_A1250_72hr.mat');
        fib_sat_limit=3;
        if ARtemp==13
            Num_points=118; % adjust number of points to ensure 4 rows of points in the discretization
        end

        %Cell geometry: general shape
        [mat_r,Npts_t,drx,dry,dA,dr_dist_squared, dist_to_line_sq, shape_name,dist_pair,drx_norm,dry_norm,Concave_ind,outline,nuc_x,nuc_y, nuc_cx, nuc_cy,choice, outside_segs, inside_segs,bdry_mat,bdry_pts,out_ind,Num_points,A,nuc_radius]=cell_geometry_units_linked_v3_multi(Num_points,A,nuc_radius,choice,nuc_decide,nuc_cx_temp,nuc_cy_temp,ARtemp);
        Tot_Integrin_Val = rho_bar*A;

        % additional parameters based on cell aspect ratio and area
        alph1=1;
        Afit=2116e-12;
        Lmax=max(max(abs(dist_pair.*inside_segs))); % maximum length scale for the cell geometry
        f_rho=(Afit/A)*103.4997/(Lmax^alph1);
        fn_tilde=(Afit/A)*0.15*(5.5e+3)/(Lmax^alph1);
        fp_tilde=(Afit/A)*0.0084*(5.5e+3)/(Lmax^alph1);

%         alph1=0;
%         Lmax=max(max(abs(dist_pair.*inside_segs))); % maximum length scale for the cell geometry
%         f_rho=1.5910e+06; %103.4997/Lmax=1.5910e+06 for square, 103.4997/Lmax=6.7561e+05 for AR=11, 103.4997/Lmax=8.4187e+05 for AR=7
%         fn_tilde= 1.2682e+07; %0.15*(5.5e+3)/Lmax=1.2682e+07 for square; 0.15*(5.5e+3)/Lmax=5.3853e+06 for AR=11; 0.15*(5.5e+3)/Lmax=6.7106e+06 for AR=7
%         fp_tilde=7.1018e+05; %0.0084*(5.5e+3)/Lmax=7.1018e+05 for square; 0.0084*(5.5e+3)/Lmax=3.0158e+05 for AR=11; 0.0084*(5.5e+3)/Lmax=3.7579e+05 for AR=7
        
        T_L=0.4;

        [L_0_rel, d_rel,EA] = bdry_tension2(dist_pair,outside_segs,mat_r,outline);
        Lambda_e=T_L.*L_0_rel+EA.*(d_rel-L_0_rel)./L_0_rel;
        % because this involves dividing by 0 in some values, some values will be
        % NaN. remove them and make them 0:
        Lambda_e(isnan(Lambda_e))=0;

        %Initial integrin concetration
        initial_choice = 1; % random free integrin initial condition
        [integrin_free, integrin_bound, rho_bar, initial_name,percent_err]=initial_integrin_units_v2_mult(Npts_t,Tot_Integrin_Val,dA,Concave_ind,outline,mat_r,rho_bar,initial_choice);

        %the maximum time: this equation determines the final time point but note
        %that this is already defined within parameter_gen_units
        time_max = time_max_factor*time_step_u;

        %initilize store variables
        time_index = 0; %time counter index
        time_points = floor(time_max/time_step_u)+1;

        % initialize the storage arrays for integrin results
        integrin_free_store = zeros(time_points,Npts_t);
        integrin_bound_store = zeros(time_points,Npts_t);
        time_store = zeros(time_points,1);
        F_store = zeros(time_points,Npts_t,2);
        %F_p_store = zeros(time_points,Npts_t,2);
        %F_n_store = zeros(time_points,Npts_t,2);
        %F_rho_store = zeros(time_points,Npts_t,2);
        R_store = zeros(time_points,Npts_t);

        % initialize the store arrays for network results
        sat_lat_idx_t = cell(1,time_points); % cell array; each element will be a vector of lattice point numbes which are saturated
        Esys_t = zeros(1,time_points); % row vect; each element will be the final Esystem value at the particular timepoint
        ControlPoints_actin_t = cell(1,time_points); % cell array; each element will be a cell containing the control points for each base fiber (combo order) at the particular timepoint
        actin_network_t = cell(1,time_points); % cell array; each element will be a cell containing the control points for each nucleus-attached fiber (combo order) at the particular timepoint
        Combo_Order_t = cell(1,time_points); % cell array; each element will be a Mx2 array identifying the combo order at the particular timepoint
        Eb_network_t = cell(1,time_points); % cell array; each element will be a 1xM array identifying the bending energies of each fiber for all combos at the particular timepoint
        Time_matrix_t = cell(1,time_points); % cell array; each element will be a matrix identifying when fibers have been constructed at each timepoint
        Lc_network_t = cell(1,time_points); % cell array; each element will be a 1xM array identifying the length of each fiber for all combos at the particular timepoint
        pt_remove_store_t = cell(1,time_points);
        NF_altered_store_t= cell(1,time_points);
        N_fibers_store_t= cell(1,time_points); % cell array: each element will be an 1xM array identifying the number of fibers that pass near any given point

        %Pp_network_t = cell(1,time_points); % cell array; each element will be a 1xM array identifying the probabilities each fiber being a premyofibril for all combos at the particular timepoint
        %Pn_network_t = cell(1,time_points); % cell array; each element will be a 1xM array identifying the probabilities each fiber being a nascent myofibril for all combos at the particular timepoint
        %Pi_p_network_t = zeros(Npts_t,Npts_t,time_points); % 3d array; each slice will be a square matrix identifying the premyofibril probability related term in the force equation each timepoint
        %Pi_n_network_t = zeros(Npts_t,Npts_t,time_points); % 3d array; each slice will be a square matrix identifying the nascent myofibril probability related term in the force equation each timepoint
        fiber_id1_t = cell(1,time_points); % cell array; each element will be a 1xM array identifying the fiber type for all combos at the particular timepoint

        % The following computations will be used repeatedly in one of the codes
        % below so we compute it here once to save computational time:
        dp_temp=dist_pair;
        dp_temp(dp_temp==0)=[]; % remove values of dist_pair that equal 0
        dA_obstruct_radius = min(dp_temp)/2; % radius of dA_obstruct = {minimum non-zero distance between 2 lattice points}/2

        % determine boundary of the cell: don't need every point on the boundary,
        % only the minimum points needed to be able to draw the boundary
        tmp=size(outline);
        if tmp(1)==2 && tmp(2)<50
            boundary_pts=outline';
        else
            % for circle, oval
            boundary_pts=bdry_pts;
        end

        StSt=0; % this term will track whether the integrin distribution has reached steady state or not. If steady state is reached, it'll be updated from 0 to 1
        StSt_thresh=0.001; % threshold value used to determine if steady state has been reached

        tmim=10; % minutes
        nt=ceil(60*tmim/time_step_u);
        lat_sat_idx=[];
        E_sys_final=0; 
        controlPoints=[]; network2=[]; combo_order=[]; Eb_network=[]; Time_matrix=[];Lc_network=[];Lat2=[];
        V2_decomp=[];pt_remove_store=[];NF_altered=[];lattice_sat_vals=[];

        for time = 0:time_step_u:time_max
            %R factors
            R_bound = R_factors_units_v2(rho_0,integrin_bound);

            % fiber network construction
            if mod(time,nt*time_step_u)==0
            [lat_sat_idx, E_sys_final, controlPoints, network2, combo_order,Eb_network,Time_matrix,Lc_network,Lat2,V2_decomp,pt_remove_store,NF_altered,lattice_sat_vals] = fiber_model_linked2_v4(mat_r,R_store,F_store,time_index,sat_lat_idx_t,...
                Combo_Order_t,kappa_limit, df, fibers_per_bund,dA_obstruct_radius,l_p,t0, boundary_pts, nuc_x,nuc_y, Fsat_lat, nuc_cx, nuc_cy, max_bundle_number,nuc_radius,time,Npts_t,Eb_max,ControlPoints_actin_t,Time_matrix_t,pt_remove_store_t,dist_pair,outline,NF_altered_store_t,StSt,dA_FA_n,fib_sat_limit,outside_segs, inside_segs,bdry_mat,choice,out_ind,nuc_rel);
            end
            % Store relevant info at this time point
            sat_lat_idx_t{time_index+1} = lat_sat_idx; % record lattice points that are currently saturated at this time point
            Esys_t(time_index+1) = E_sys_final; % record final Esystem value 
            ControlPoints_actin_t{time_index+1} = controlPoints; % record control points used to create base fibers at this time point
            actin_network_t{time_index+1} = network2; % record control points used to create nucleus-attached fibers at this time point
            Combo_Order_t{time_index+1} = combo_order; % record combo order at this time point
            Eb_network_t{time_index+1} = Eb_network;
            Time_matrix_t{time_index+1} =Time_matrix;
            Lc_network_t{time_index+1} = Lc_network;
            pt_remove_store_t{time_index+1}=pt_remove_store;
            NF_altered_store_t{time_index+1}=NF_altered;
            N_fibers_store_t{time_index+1}=lattice_sat_vals;

            % Determine fiber identification and associated probabilities for actin
            % network
            % Grab force vectors at each point: this is basically a Npts_tx2 array
            F_vectors = F_store(max(time_index-1,1),:,:);  % 
            F_x = F_vectors(1,:,1); %F_x is the 1st column of the force vector; doing this creates a ROW vector
            F_y = F_vectors(1,:,2);%F_y is the 2nd column of the force vector; doing this creates a ROW vector    

            [Pp_network, Pn_network,fiber_id1] = fiber_id_linked2(Time_matrix_t,Npts_t, Eb_max, fibers_per_bund,time_index,a,b,F_x, F_y, Fsat_lat,F_c);
            %Pp_network_t{time_index+1}=Pp_network;
            %Pn_network_t{time_index+1}=Pn_network;
            fiber_id1_t{time_index+1}=fiber_id1;

            %Force
            [F, F_adh, F_cyto, F_p, F_n] = Force_units_Lavg_linked_v3(R_bound,Npts_t,f_rho,drx,dry,dA,rho_sat,Lc_network,Time_matrix,lat_sat_idx,V2_decomp,Pp_network,Pn_network,fp_tilde,fn_tilde,dA_FA_p,dA_FA_n,rho_0, outside_segs, inside_segs,bdry_mat,drx_norm,dry_norm,Lambda_e);
            %Pi_p_network_t(:,:,time_index+1)=Pi_p_actin_network;
            %Pi_n_network_t(:,:,time_index+1)=Pi_n_actin_network;

            %store values for this time step - note that this is the state of the
            %previous time step -- this means that I don't store the values for the
            %last time step, which shouldn't produce much of a problem
            time_index = time_index+1;
            integrin_free_store(time_index,:) = integrin_free;
            integrin_bound_store(time_index,:) = integrin_bound;
            time_store(time_index) = time;
            F_store(time_index,:,:) = F;
            %F_p_store(time_index,:,:) = Fp;
            %F_n_store(time_index,:,:) = Fn;
            R_store(time_index,:) = R_bound;

            %next time step integrin cocentrations    
            [integrin_free_new, integrin_bound_new]=New_integrin_cons_units_v2(integrin_free, integrin_bound,F, k_0,k_1,k_m1,time_step_u,rho_bar,Concave_ind,F_0);

            % code error checkpoint
            if max(isnan(integrin_free_new)) >0
                disp('ERROR! Free integrin density is unbounded')
                break;
            end
            %It is possible for the code to vear off away from a proper solution,
            %if that is the case the integrin concentration becomes non-real. The
            %code will stop if that is the case.
            if max(isnan(integrin_bound_new)) >0
                disp('ERROR! Bound integrin density is unbounded')
                break;
            end
            if any(integrin_bound_new<0) 
                disp('ERROR! Bound integrin is predicted to be negative')
                break;
            end    
            integrin_free_old=integrin_free;
            integrin_bound_old=integrin_bound;

            %If the time step is chosen to be too large, it is possible to loose
            %mass consrvation. The following is set up to stop the simulation is
            %that is the case
            if (dA*sum(integrin_free_new+integrin_bound_new)...
                    < 0.999*Tot_Integrin_Val)||(dA*sum(integrin_free_new+integrin_bound_new) > 1.001*Tot_Integrin_Val)
                disp('ERROR there is no mass conservation exiting code. \\ Think of decreasing step size')
                disp(num2str(time))
                break;
            end
            %reset previous time step to new time step
            integrin_bound=integrin_bound_new;
            integrin_free=integrin_free_new;

            %to keep track of how the code is running we display every n steps
            if rem(time_index,30)==0
                disp([num2str(time) ' time step out of ' num2str(time_max)])
            end

            % Determine if steady state has been reached in the integrin
            % distibution: if the percent change is less than StSt_thresh, then
            % steady state is reached
            integrin_tot_old=integrin_free_old+integrin_bound_old;
            integrin_tot_new=integrin_free_new+integrin_bound_new;
            integrin_diff=norm((integrin_tot_new-integrin_tot_old)./integrin_tot_old);
            if integrin_diff <= StSt_thresh
                StSt=1;
            end
        end
        
        % store network results at each hour (storing for every time point may
        % result in an extremely large file so we'll store only the network at each
        % hour)
        thr=1:floor(t_max_u/3600);
        fiber_id1_t_store=cell(1,numel(thr));
        actin_network_t_store=cell(1,numel(thr));
        Combo_Order_t_store=cell(1,numel(thr));
        tsec=thr*3600; % convert from hours to seconds 
        tp=floor(1+tsec/time_step_u);
        for j=1:numel(thr)
            fiber_id1_t_store{j}=fiber_id1_t{tp(j)};
            actin_network_t_store{j}=actin_network_t{tp(j)};
            Combo_Order_t_store{j}=Combo_Order_t{tp(j)};
        end

        % store network results at each hour
        net_res{my2,my1}=actin_network_t_store; %actin_network_t{end};
        F_res{my2,my1}=F_store;
        fiber_id_res{my2,my1}=fiber_id1_t_store; %fiber_id1_t{end};
        combo_res{my2,my1}=Combo_Order_t_store; %Combo_Order_t{end};
        nuc_cx_store(my2,my1)=nuc_cx;
        nuc_cy_store(my2,my1)=nuc_cy;
        outline_store{my2,my1}=outline;
    end
    dA_store(my2)=dA;
end
    

%save('square_fsat_bet_param_ar=9-2_Ntot_results_store','net_res','combo_res','F_res','fiber_id_res','mat_r','nuc_cx_store','nuc_cy_store','dA');

filename_save = ['paramTEST5_A1250_72hr_fsatlim=3_net2_results_store.mat'];
save(filename_save,'net_res','combo_res','F_res','fiber_id_res','mat_r','nuc_cx_store','nuc_cy_store','dA_store','outline_store','nuc_radius','nsim','AR_vec');

%Display the time it took to run the code in hours and minutes
elapsed_time_hrs = floor(toc/(60.*60.));
elapsed_time_min = floor((toc - 60*60*elapsed_time_hrs)/60);
elapsed_time_sec = toc - 60*elapsed_time_min- 60*60*elapsed_time_hrs;
disp(['It took ' num2str(elapsed_time_hrs) ' hrs and ' num2str(elapsed_time_min) ' min to run the code']);
