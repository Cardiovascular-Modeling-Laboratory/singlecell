% this code runs single_cell_units_linked_v3 multiple times on a square
% cell with the nucleus placement altered.
% in the resuls, each row corresonds to one of the nucleus placements:

clear;

% Define material properties (outside the loop)
E = 1000; % Young's modulus in Pascals (Pa)
Poisson_ratio = 0.5;  % Poisson's ratio
% Define stretch schedule as a matrix [start_time, end_time, stretch_factor_start, stretch_factor_end]
stretch_schedule = [
    0,     24*3600, 1.0, 1.0;   % 0 h - 24 h: No stretch
    24*3600, 48*3600, 1.0, 1.5; % 24 h - 48 h: Stretch from 1x to 5x
    48*3600, 72*3600, 1.5, 1.0; % 48 h - 72 h: Shrink back to 1x
    ];
nsim=1; % max number of simulations
choice=1; % designate the square
nuc_decide_vec=[1 1 3 3]; % vector indicating if the nucleus should be placed in the center (1), randomly (2), or manually (3)
nuc_cx_temp_vec=[0 0 -7 0]; % vector indicating how x-coord of nucleus should be offset if being placed manually
nuc_cy_temp_vec=[0 0 0 7]; % vector indicating how y-coord of nucleus should be offset if being placed manually
nuc_rel_vec=[0 1 1 1]; %[0 1 1 1]; % vector indicating if nucleus is treated as an obstruction (1) or not (0)
fib_sat_limit_vec=[0 18 18 18]; % vector indicating base saturation number of lattice points within the nucleus
%fib_sat_limit_vec=[0 3 3 3]; % vector indicating base saturation number of lattice points within the nucleus

net_res=cell(numel(nuc_rel_vec),nsim);
F_res=cell(numel(nuc_rel_vec),nsim);
fiber_id_res=cell(numel(nuc_rel_vec),nsim);
combo_res=cell(numel(nuc_rel_vec),nsim);
nuc_cx_store=zeros(numel(nuc_rel_vec),nsim);
nuc_cy_store=zeros(numel(nuc_rel_vec),nsim);

tic



for my2=1:numel(nuc_rel_vec)
    disp(['nucleus placement number ',num2str(my2)])
    nuc_rel=nuc_rel_vec(my2);
    nuc_decide=nuc_decide_vec(my2);
    nuc_cx_temp=nuc_cx_temp_vec(my2);
    nuc_cy_temp=nuc_cy_temp_vec(my2);
    
    for my1=1:nsim
        disp(['  sim number ',num2str(my1)])
        %load parameters
        load('paramTEST5_72hr.mat');
        fib_sat_limit=fib_sat_limit_vec(my2);
        
        %Cell geometry: general shape
        % TODO: Currently, we initialize the cell shape only once at the beginning of the loop
        % What we want to do is update the cell shape at each time step to simulate stretching
        [mat_r,Npts_t,drx,dry,dA,dr_dist_squared, dist_to_line_sq, shape_name,dist_pair,drx_norm,dry_norm,Concave_ind,outline,nuc_x,nuc_y, nuc_cx, nuc_cy,choice, outside_segs, inside_segs,bdry_mat,bdry_pts,out_ind,Num_points,A,nuc_radius]=cell_geometry_units_linked_v3_nuc(Num_points,A,nuc_radius,choice,nuc_decide,nuc_cx_temp,nuc_cy_temp);
        
        disp(['mat_r: ']);
        mat_r
        % additional parameters based on cell aspect ratio and area
        alph1=1;
        Lmax=max(max(abs(dist_pair.*inside_segs))); % maximum length scale for the cell geometry
        f_rho=103.4997/(Lmax^alph1);%134.7428/(Lmax^alph1);
        fn_tilde=0.15*(5.5e+3)/(Lmax^alph1);%0.154*(5.5e+3)/(Lmax^alph1);%0.2*(5.5e+3)/(Lmax^alph1);
        fp_tilde=0.0084*(5.5e+3)/(Lmax^alph1);%0.0111*(5.5e+3)/(Lmax^alph1);
        
        T_L=0.4;
        
        [L_0_rel, d_rel,EA] = bdry_tension2(dist_pair,outside_segs,mat_r,outline);
        Lambda_e=T_L.*L_0_rel+EA.*(d_rel-L_0_rel)./L_0_rel;
        % because this involves dividing by 0 in some values, some values will be
        % NaN. remove them and make them 0:
        Lambda_e(isnan(Lambda_e))=0;
        
        %Initial integrin concetration
        initial_choice = 2; % random free integrin initial condition
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
        
        % Parameters for cell geometry modifications
        initial_area = dA * Npts_t;
        initial_dA = dA;
        A_initial = A;
        dA_old = dA;
        Npts_t_old = Npts_t;
        mat_r_initial = mat_r;
        outline_initial = outline;
        nuc_cx_initial = nuc_cx;
        nuc_x_initial = nuc_x;
        nuc_cy_initial = nuc_cy;
        nuc_y_initial = nuc_y;
        % Initialize arrays to store results for later analysis
        F_t_store = zeros(time_points, 1);                % Expected force (N)
        F_adhesion_store = zeros(time_points, 1);         % Adhesion force from simulation (N)
        force_difference_store = zeros(time_points, 1);   % Difference between adhesion and expected force (N)
        force_ratio_store = zeros(time_points, 1);        % Ratio of adhesion force to expected force (dimensionless)
        
        % Calculate initial length L_0 (before the loop)
        x_positions = mat_r_initial(:,1);
        L_0 = (max(x_positions) - min(x_positions)) * 1e-6; % Convert to meters if necessary
        
        for time = 0:time_step_u:time_max
            %R factors
            R_bound = R_factors_units_v2(rho_0,integrin_bound);
            
            % Simulate stretching of the cell
            % Get the stretch factor based on the schedule and current time
            stretch_factor = get_stretch_factor(time, stretch_schedule);
            disp(['Time: ', num2str(time), ' Stretch factor: ', num2str(stretch_factor)]);
            % Recalculate cell geometry
            % If loop number is > 1
            if time >= 0
                
                [mat_r,Npts_t,drx,dry,dA,dr_dist_squared, dist_to_line_sq, shape_name,dist_pair,drx_norm,dry_norm,Concave_ind,outline,nuc_x,nuc_y,nuc_cx,nuc_cy,choice, outside_segs, inside_segs,bdry_mat,bdry_pts,out_ind,Num_points,A,nuc_radius] = cell_geometry_strecth_v1(A_initial, stretch_factor, mat_r_initial, Npts_t,drx,dry,dA,dr_dist_squared, dist_to_line_sq, shape_name,dist_pair,drx_norm,dry_norm,Concave_ind,outline_initial,nuc_x_initial,nuc_y_initial,nuc_cx_initial,nuc_cy_initial,choice, outside_segs, inside_segs,bdry_mat,bdry_pts,out_ind,Num_points,A,nuc_radius,Poisson_ratio);
                tmp=size(outline);
                if tmp(1)==2 && tmp(2)<50
                    boundary_pts=outline';
                else
                    % for circle, oval
                    boundary_pts=bdry_pts;
                end
                
                % % Recalculate Lambda_e
                [L_0_rel, d_rel, EA] = bdry_tension2(dist_pair, outside_segs, mat_r, outline);
                Lambda_e = T_L .* L_0_rel + EA .* (d_rel - L_0_rel) ./ L_0_rel;
                Lambda_e(isnan(Lambda_e)) = 0;
                
                % インテグリン濃度の調整
                disp(['Adjusting integrin concentration for cell geometry stretching at time ', num2str(time)]);
                disp(['Stretch factor: ', num2str(stretch_factor)]);
                % disp([ 'Integrin free: ', num2str(sum(integrin_free)), ' Integrin bound: ', num2str(sum(integrin_bound))]);
                % disp(['dA = ', num2str(dA), ' Npts_t = ', num2str(Npts_t)]);
                area_new = dA * Npts_t;
                area_scale_factor = initial_area / area_new;
                % integrin_free = integrin_free * area_scale_factor;
                % integrin_bound = integrin_bound * area_scale_factor;
                % disp(['area_new = ', num2str(area_new), ' initial_area = ', num2str(initial_area)]);
                % disp(['Area scale factor: ', num2str(area_scale_factor)]);
                % disp([ 'Integrin free: ', num2str(sum(integrin_free)), ' Integrin bound: ', num2str(sum(integrin_bound))]);
                % disp(['dA*sum(integrin_free_new+integrin_bound_new) =', num2str(dA*sum(integrin_free+integrin_bound))]);
                
                % 次のイテレーションのために古い値を更新
                dA_old = dA;
                Npts_t_old = Npts_t;
            end
            
            % Check if the cell geometry is correctly updated
            disp(['Time: ', num2str(time), ' Cell geometry updated']);
            disp(['dA: ', num2str(dA), ' Npts_t: ', num2str(Npts_t), ' Area: ', num2str(dA * Npts_t)]);
            % mat_r
            disp(['drx: ', num2str(size(drx)), ' dry: ', num2str(size(dry))]);
            
            % fiber network construction
            if mod(time,nt*time_step_u)==0
                disp(['***']);
                disp(['Time: ', num2str(time), ' Fiber network construction']);
                disp(['***']);
                [lat_sat_idx, E_sys_final, controlPoints, network2, combo_order,Eb_network,Time_matrix,Lc_network,Lat2,V2_decomp,pt_remove_store,NF_altered,lattice_sat_vals] = fiber_model_linked2_v4(mat_r,R_store,F_store,time_index,sat_lat_idx_t,...
                    Combo_Order_t,kappa_limit, df, fibers_per_bund,dA_obstruct_radius,l_p,t0, boundary_pts, nuc_x,nuc_y, Fsat_lat, nuc_cx, nuc_cy, max_bundle_number,nuc_radius,time,Npts_t,Eb_max,ControlPoints_actin_t,Time_matrix_t,pt_remove_store_t,dist_pair,outline,NF_altered_store_t,StSt,dA_FA_n,fib_sat_limit,outside_segs, inside_segs,bdry_mat,choice,out_ind,nuc_rel);
            end
            
            % Check if actin_network is populated
            disp(['Time: ', num2str(time), ' Actin network size: ', num2str(size(network2))]);
            disp(['Actin network: ']);
            network2
            disp(['Esys']);
            E_sys_final
            disp(['Control points']);
            controlPoints
            disp(['Combo order']);
            combo_order
            
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
            
            
            % After updating F_store and time_index
            [F_t, F_adhesion, force_difference, force_ratio] = ...
                geometry_change_force_validation_v1(E, Poisson_ratio, stretch_factor, initial_area, L_0, ...
                F_store, time_index, dA, time);
            % Optionally, store the results for later analysis
            F_t_store(time_index) = F_t;
            F_adhesion_store(time_index) = F_adhesion;
            force_difference_store(time_index) = force_difference;
            force_ratio_store(time_index) = force_ratio;
            
            
            %next time step integrin cocentrations
            % Print the amount of integrin before and after the time step
            % disp(['**']);
            % disp(['Integrin free: ', num2str(sum(integrin_free)), ' Integrin bound: ', num2str(sum(integrin_bound)), ' rho_bar: ', num2str(rho_bar)]);
            [integrin_free_new, integrin_bound_new]=New_integrin_cons_units_v2(integrin_free, integrin_bound,F, k_0,k_1,k_m1,time_step_u,rho_bar,Concave_ind,F_0);
            % For now, we apply area_scale_factor to integrin concentrations
            % if area_scale_factor is defined (i.e. the cell shape has changed)
            if exist('area_scale_factor', 'var')
                % integrin_free_new=integrin_free_new * area_scale_factor;
                % integrin_bound_new=integrin_bound_new * area_scale_factor;
            end
            % disp(['Integrin free new: ', num2str(sum(integrin_free_new)), ' Integrin bound new: ', num2str(sum(integrin_bound_new)), ' rho_bar: ', num2str(rho_bar)]);
            % disp(['**']);
            
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
            % TODO: Modify this to allow cell geometry streching
            disp(['integrin_free: ', num2str(sum(integrin_free)), ' integrin_bound: ', num2str(sum(integrin_bound))]);
            disp(['integrin_free_new+: ', num2str(sum(integrin_free_new)), ' integrin_bound_new+: ', num2str(sum(integrin_bound_new))]);
            disp(['integrin_free + integrin_bound: ', num2str(sum(integrin_free+integrin_bound))]);
            disp(['dA*sum(integrin_free_new+integrin_bound_new): ', num2str(dA*sum(integrin_free_new+integrin_bound_new))]);
            disp(['Total Integrin Val: ', num2str(Tot_Integrin_Val)]);
            % disp(['dA = ', num2str(dA), ' Npts_t = ', num2str(Npts_t), ' Total Integrin Val: ', num2str(Tot_Integrin_Val)]);
            % disp(['dA*sum(integrin_free_new+integrin_bound_new) = ', num2str(dA*sum(integrin_free_new+integrin_bound_new))]);
            % disp(['0.999*Tot_Integrin_Val = ', num2str(0.999*Tot_Integrin_Val)]);
            % disp(['1.001*Tot_Integrin_Val = ', num2str(1.001*Tot_Integrin_Val)]);
            % disp(['************']);
            if ~exist("area_scale_factor", "var")
                if (dA*sum(integrin_free_new+integrin_bound_new)...
                        < 0.999*Tot_Integrin_Val)||(dA*sum(integrin_free_new+integrin_bound_new) > 1.001*Tot_Integrin_Val)
                    disp('ERROR there is no mass conservation exiting code. \\ Think of decreasing step size')
                    disp(num2str(time))
                    break;
                end
            end
            %reset previous time step to new time step
            integrin_bound=integrin_bound_new;
            integrin_free=integrin_free_new;
            
            %to keep track of how the code is running we display every n steps
            if rem(time_index,30)==0
                disp([num2str(time) ' time step out of ' num2str(time_max), ' Stretch Factor: ', num2str(stretch_factor)])
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
            
            % TODO: change the cell geometry here to integrate new feature: Dynamically change the cell shape to imitate stretching
        end
        
        % store network results at each hour (storing for every time point may
        % result in an extremely large file so we'll store only the network at each
        % hour)
        % BUG: Results not stored
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
        
        % For force validation
        F_t_res{my2,my1}=F_t_store;
        F_adhesion_res{my2,my1}=F_adhesion_store;
        force_difference_res{my2,my1}=force_difference_store;
        force_ratio_res{my2,my1}=force_ratio_store;
    end
end

%base_filename = file(1:(length(file)-4));
% Add datetime_ms to the filename
filename_save = ['paramTEST5_72hr_NucPlacement_fsatlim=18_net2_results_store_', datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.mat'];
save(filename_save,'net_res','combo_res','F_res','fiber_id_res','mat_r','nuc_cx_store','nuc_cy_store','dA','outline','nuc_radius','nsim','F_t_res','F_adhesion_res','force_difference_res','force_ratio_res','time_store','stretch_schedule');


%Display the time it took to run the code in hours and minutes
elapsed_time_hrs = floor(toc/(60.*60.));
elapsed_time_min = floor((toc - 60*60*elapsed_time_hrs)/60);
elapsed_time_sec = toc - 60*elapsed_time_min- 60*60*elapsed_time_hrs;
disp(['It took ' num2str(elapsed_time_hrs) ' hrs and ' num2str(elapsed_time_min) ' min to run the code']);


% Function to compute stretch factor based on the current time and schedule
function stretch_factor = get_stretch_factor(time, stretch_schedule)
    for i = 1:size(stretch_schedule, 1)
        start_time = stretch_schedule(i, 1);
        end_time = stretch_schedule(i, 2);
        factor_start = stretch_schedule(i, 3);
        factor_end = stretch_schedule(i, 4);
        
        if time >= start_time && time <= end_time
            % Linearly interpolate stretch factor within the time interval
            stretch_factor = factor_start + (factor_end - factor_start) * ((time - start_time) / (end_time - start_time));
            return;
        end
    end
    % Default stretch factor if no match (e.g., outside the schedule)
    stretch_factor = 1.0;
end
