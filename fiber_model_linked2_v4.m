function [lat_sat_idx, E_sys_final, controlPoints, network2, combo_order,Eb_network,Time_matrix,Lc_network,Lat2,V2_decomp,pt_remove_store,NF_altered,lattice_sat_vals] = fiber_model_linked2_v4(mat_r,R_store,F_store,time_index,sat_lat_idx_t,...
    Combo_Order_t,kappa_limit, df, fibers_per_bund,dA_obstruct_radius,l_p,t0, boundary_pts, nuc_x,nuc_y, Fsat_lat, nuc_cx, nuc_cy, max_bundle_number,nuc_radius,time,Npts_t,Eb_max,ControlPoints_actin_t,Time_matrix_t,pt_remove_store_t,dist_pair,outline,NF_altered_store_t,StSt,dA_FA_n,fib_sat_limit,outside_segs, inside_segs,bdry_mat,choice,out_ind,nuc_rel)

tsec_idx=time_index+1; % use time_index+1 since time_index starts at 0

% This code was written using the micometer scale to simplify computations.
% So, convert all relevant parameters/terms from meter to micrometer
outline=10^6.*outline;
Lat=10^6.*mat_r; % mat_r is in units of meter, multiply by 10^6 to put it in units of micrometer
dA_obstruct_radius=10^6.*dA_obstruct_radius; % convert from meters to micrometers
nuc_radius=10^6*nuc_radius;
nuc_cx=10^6*nuc_cx;
nuc_cy=10^6*nuc_cy;
nuc_x=10^6.*nuc_x;
nuc_y=10^6.*nuc_y;
boundary_pts=10^6.*boundary_pts;
dA_FA_n=dA_FA_n*10^12; 

%pt_remove=unique([0 pt_remove_store_t{max(tsec_idx-1,1)}]);
pt_remove=0;
NF_altered=[0 0];
%NF_altered=unique([0 0; NF_altered_store_t{max(tsec_idx-1,1)}],'rows');
% if integrin steady state has been reached, then we can cut down on
% computational time by providing the information needed to build the
% steady state network
% if StSt==1
%     pt_remove=unique([0 pt_remove_store_t{max(tsec_idx-1,1)}]);
%     NF_altered=NF_altered_store_t{max(tsec_idx-1,1)};
% end

ncnts=1;
pt_remove_store=pt_remove;

while ~isempty(pt_remove) && ncnts<=1000 && ~isempty(NF_altered)
    % arrays for storage
    Eb_network = cell(Npts_t,Npts_t); % array to track/store bending energy of each constructed fiber
    Lc_network = cell(Npts_t,Npts_t); % array to track/store bending energy of each constructed fiber
    Time_matrix=[]; % array to track/store the info needed to construct the probability function of each constructed fiber (used in a later code)
    lattice_sat_vals_temp=[]; % initialize this to be empty; if no points are saturated, this will allow the nesprin portion to still run
    lattice_sat_vals=zeros(size(Lat(:,1)));
    
    % Grab force vectors at each point: this is basically a Npts_tx2 array
    F_vectors = F_store(max(tsec_idx-1,1),:,:);  % 
    F_x = F_vectors(1,:,1); %F_x is the 1st column of the force vector; doing this creates a ROW vector
    F_y = F_vectors(1,:,2);%F_y is the 2nd column of the force vector; doing this creates a ROW vector
    % F_x and F_y are row vectors but we want the force array to have 2 columns
    % where 1st column is Fx and 2nd column is Fy. To do this, create an array
    % that has 2 rows (Fx as row 1, Fy as row 2) and then just transpose the
    % array:
    V1=[F_x; F_y];
    V=V1'; % this is the array where the 1st (2nd) column is the x (y) component of the force vector at each lattice point

    % Isolate only the lattice points which are relevant
    [lat_sat_idx,lat_nuc_idx,~] = isolatedLattice3_linked_v2(Lat,nuc_radius,nuc_cx,nuc_cy,Fsat_lat,V,nuc_rel,R_store,tsec_idx) ;
    lattice_sat_vals(lat_nuc_idx)=fib_sat_limit;
    % there must be at least 2 lattice saturated lattice points in order
    % for the network construction to begin. So check that this is the case
    if numel(lat_sat_idx)<2
        E_sys_final=0;
        controlPoints=cell(1,1);
        network2=cell(1,1);
        combo_order=[];
        Lat2 = Lat(lat_sat_idx,:);
        V2_decomp=cell(1,1);
        pt_remove_store=[];
        return
    end
    lat_sat_idx(ismember(lat_sat_idx,pt_remove_store))=[];
    Lat2 = Lat(lat_sat_idx,:); % lattice of points which are saturated

    % isolate the corresponding tangent vectors for the starting/ending points
    V2 = V(lat_sat_idx,:);
    R_vals=R_store(max(1,tsec_idx-1),lat_sat_idx);
    c=1.233;
    vec_decomp_numbs=2.*floor(((R_vals+0.05)./(c-1)).^(1/c))-1; % number of vectors to use when decomposing each force vector

    % only decompose vectors which are sufficiently angled. If they're not
    % sufficiently angled, they are not decomposed so set their decomp
    % value to 1
    ad=atand(V2(:,2)./V2(:,1));
    deg_max=8;%15; 
    vec_decomp_numbs(find(abs(ad)<=deg_max))=1;
    
    % check if the number of vectors used at any point needs to be altered
    altered_pts_idx=intersect(lat_sat_idx,NF_altered(:,1));
    if ~isempty(altered_pts_idx)
        % identify the points
        reduced_idx=find(ismember(lat_sat_idx,altered_pts_idx));
        reduced_idx2=find(ismember(NF_altered(:,1),altered_pts_idx));

        % update vec_decomp_numbs at the designated pts
        vec_decomp_numbs(reduced_idx)=NF_altered(reduced_idx2,2);
    end

    % determine number of vectors each lattice point can be decomposed into
    % and record the vector decomposition
    V2_decomp=cell(1,size(V2,1));
    vec_track=zeros(size(V2,1),1); % track the number of decomposition vectors 
    N_F_max=5; % maximum number of vectors allowed in the force decomposition
    for vcomp=1:size(V2,1)
        v=V2(vcomp,:);
        pt=Lat2(vcomp,:);
        N_F=vec_decomp_numbs(vcomp);
        [W, W_n] = v_splitting3(v,N_F,N_F_max,pt,outline,Lat,mat_r,choice,out_ind,outside_segs,lat_sat_idx); 
        V2_decomp{vcomp}=W;
        vec_track(vcomp)=W_n;
    end
    V2_unused=V2_decomp; % this is a cell array containing all the unused vectors; initially, it is just the array of all possible vectors

    numberOfPts=size(Lat2,1); % number of saturated lattice points

        if numel(lat_sat_idx)>1 
            [~, A_template, combo_order]=network_template(Lat,outline, V2, Lat2,outside_segs, inside_segs,bdry_mat,lat_sat_idx,choice,vec_decomp_numbs);
        elseif numel(lat_sat_idx)==1 
            combo_order=[lat_sat_idx, lat_sat_idx];
        else
            combo_order=[];
        end
%        numbOfCombs=size(combo_order,1);
%     % Update/Set combo order
%     if tsec_idx==1 % if this is the first simulation, construct pt-pt combinations 
%         % construct the pt-pt combinations
%         if numel(lat_sat_idx)>1 
%             [~, A_template, combo_order]=network_template(Lat,outline, V2, Lat2,outside_segs, inside_segs,bdry_mat,lat_sat_idx);
%         elseif numel(lat_sat_idx)==1 
%             combo_order=[lat_sat_idx, lat_sat_idx];
%         else
%             combo_order=[];
%         end
%         numbOfCombs=size(combo_order,1);        
%     else % else, pt-pt combos have already been constructed at a previous time point. So just adjust the combo order based on how the saturated lattice has changed
%         sat_lattice_old=sat_lat_idx_t{tsec_idx-1};
%         sat_lattice_new=lat_sat_idx;
%         % determine which subcase we're in: is the lattice the same or have
%         % points been added/removed?
%         case1=isequal(sat_lattice_old,sat_lattice_new); % =1 if lattices are identical, =0 if not identical
%         switch case1
%             case 1 % if the 2 lattices are identical, keep combo order the same as before
%                 [~, A_template, combo_order]=network_template(Lat,outline, V2, Lat2,outside_segs, inside_segs,bdry_mat,lat_sat_idx,choice,vec_decomp_numbs);
%                 % unless we're constructing from nothing, maintain combo
%                 % order
%                 if ~isempty(sat_lattice_new)
%                     combo_order_old = Combo_Order_t{tsec_idx-1};
%                     combo_order = combo_order_old;
%                 end
%             case 0 % the lattices are not identical. So keep the points which have been maintained. All other points should be randomized and appended to the "maintained points" order list
%                 [~, A_template, combo_order]=network_template(Lat,outline, V2, Lat2,outside_segs, inside_segs,bdry_mat,lat_sat_idx,choice,vec_decomp_numbs);
%         end
        if ~isempty(combo_order)
            idx_remove=union(find(ismember(combo_order(:,1),pt_remove_store)==1),find(ismember(combo_order(:,2),pt_remove_store)==1));
            combo_order(idx_remove,:)=[];
        end
        numbOfCombs=size(combo_order,1);
%    end

    % create storage for results
    network=cell(numbOfCombs,1); % cell array that will hold the r_x(t) and r_y(t) for each bundle
    E_bundle = cell(numbOfCombs,1); % total bundle bending energy for each pt-pt combo
    controlPoints=cell(numbOfCombs,1);
    E_sys_tot = zeros(numbOfCombs,1);

    %disp(['time: ',num2str(time),'; There are ', num2str(numberOfPts), ' saturated lattice points and ', num2str(numbOfCombs),' potential point-point combinations'])
    %disp(['time: ',num2str(time), '; There are ', num2str(numberOfPts), ' saturated lattice points '])

    Currentopenpts=ones(numberOfPts,1); % a 1 indiciates that that lattice point is open; a 0 will indicate that it is occupied

    % Main code: curve construction
    % keep track of whether vectors are being used in each combo. 1
    % indicates it has been used to construct a fiber; 0 indicates the
    % vector was not used to construct a fiber
    vecs_used=cell(1,numberOfPts);
    for vind=1:numberOfPts
        vecs_used{vind}=zeros(1,vec_track(vind));
    end
    t=linspace(0,1);
    
    if numbOfCombs>0
        [costs,~] = dijkstra(A_template,Lat2); % energetic cost associated with the template network
    end

    for combo_number=1:size(combo_order,1)
        starting_pt=combo_order(combo_number,1);
        ending_pt=combo_order(combo_number,2);
        P0 = Lat(starting_pt,:); % record starting lattice point
        P4 = Lat(ending_pt,:); % record ending lattice point

        [proceedq, int_idx]= proceed_id(tsec_idx,Combo_Order_t,starting_pt,ending_pt); 
        pt1=find(lat_sat_idx==starting_pt);
        pt2=find(lat_sat_idx==ending_pt);
        W1_temp=V2_decomp{pt1}; % array of decomposition vectors for point 1
        W2_temp=V2_decomp{pt2}; % array of decomposition vectors for point 2
        rt_vec_storage=cell(1,size(W1_temp,1));
        Eb_temp_store=-1.*ones(N_F_max,max_bundle_number); % temporary storage matrix for bending energies
        Lc_temp_store=-1.*ones(N_F_max,max_bundle_number); % temporary storage matrix for curve lengths

        Time_info_new=[];
        controlPoints{combo_number}=cell(N_F_max,max_bundle_number);

        for vec_ctr=1:size(W1_temp,1) %usually size(W1)temp,1)=6 but it may =1 if vector is horizontal or vertical, i.e., the number of decomposition vectors I've imposed
            % attempt to place curve bundles
            bundle_counter=1;  % bundle number counter
            keep_curve=1; % initialize to enter while loop

            rt=zeros(2,length(t),max_bundle_number); % this is a 3 dimensional array which 
            % will hold the r(t) information for each fiber.
            % It has 2 rows (for x(t) and y(t)), length(t) columns (because there are
            % length(t) coordinate points) and this is done combNumber of times (r(t)
            % is computed for each combination)        

            W1=W1_temp(vec_ctr,:);

            while keep_curve>0 && bundle_counter<= max_bundle_number && ~isempty(W2_temp)
                % Given P0,P4 we can find the best pair of subvectors w1i and
                % w2j from the tangent vector decompositions    
                W2=W2_temp;
                [w1i,w2j,I_col]=wiwjbest3(P0,P4,W1,W2);

                % CURVE CONSTRUCTION
                if bundle_counter==1 % construct initial curve 
                    % construct min energy curve, assuming curve avoids obstructions
                    [P1,P2,P3,Fiber_Bending_Energy_std,rtx,rty] = MinEnergyCurveWithoutObstruct3(P0,P4,w1i,w2j,t,kappa_limit,proceedq,ControlPoints_actin_t,tsec_idx,int_idx,vec_ctr);

                    % check that this initial proposed curve avoids the nucleus. If
                    % it doesn't, construct one that does
                    hits_nuc = avoid_nuc(rtx,rty,nuc_x,nuc_y,nuc_rel);
                    if hits_nuc==1 % note: this will produce 2 possible curves
                        [P2_1,P2_2,Fiber_Bending_Energy_std_1,Fiber_Bending_Energy_std_2,rtx_1,rty_1,rtx_2,rty_2] = MinEnergyCurveWithObstruct(P0,P1,P3,P4,t,nuc_cx,nuc_cy,nuc_radius);
                    else
                        P2_1=P2;
                        P2_2=P2;
                        Fiber_Bending_Energy_std_1=Fiber_Bending_Energy_std;
                        Fiber_Bending_Energy_std_2=Fiber_Bending_Energy_std;
                        rtx_1=rtx;
                        rty_1=rty;
                        rtx_2=rtx;
                        rty_2=rty;
                    end
                    % if a curve for this pt-pt combo has been placed, then attempt to create more curves via nudging
                else % construct nudged curve
                    P1=controlPoints{combo_number}{vec_ctr,bundle_counter-1}(2,:);
                    P2=controlPoints{combo_number}{vec_ctr,bundle_counter-1}(3,:);
                    P3=controlPoints{combo_number}{vec_ctr,bundle_counter-1}(4,:);
                    [P1,P2,P3,Fiber_Bending_Energy_std,rtx,rty] = nudged_curve(P0,P1,P2,P3,P4,1.8*t0,t); 
                    
                    P2_1=P2;
                    P2_2=P2;
                    Fiber_Bending_Energy_std_1=Fiber_Bending_Energy_std;
                    Fiber_Bending_Energy_std_2=Fiber_Bending_Energy_std;
                    rtx_1=rtx;
                    rty_1=rty;
                    rtx_2=rtx;
                    rty_2=rty;
                end % end if statement for curve construction for the P0-P4 combo
                
                % create temporary storage to identify best curve for
                % placement
                Rtx=[rtx; rtx_1; rtx_2];
                Rty=[rty; rty_1; rty_2];
                P2_temp=[P2; P2_1; P2_2];
                Fiber_Bending_Energy_std_temp=[Fiber_Bending_Energy_std; Fiber_Bending_Energy_std_1; Fiber_Bending_Energy_std_2];
                
                E_current_sys_temp=E_sys_tot(combo_number).*ones(3,1); % energy of the current system
                E_prop_sys=zeros(3,1);
                E_new_sys=zeros(3,1);
                place_curve_store=zeros(3,1);
                curve_prob_store=zeros(3,1);
                keep_curve_store=zeros(3,1);
                FiberBendingEnergy_curve_store=zeros(3,1);
                lattice_sat_vals_test=cell(3,1);
                
                for esys_idx=1:numel(Fiber_Bending_Energy_std_temp)
                    rtx_test=Rtx(esys_idx,:);
                    rty_test=Rty(esys_idx,:);
                    P2_test=P2_temp(esys_idx,:);
                    Fiber_Bending_Energy_std_test=Fiber_Bending_Energy_std_temp(esys_idx,:);
                    
                    % CURVE PLACEMENT
                    [FiberBendingEnergy_curve, keep_curve,curve_prob]=curve_comps(rtx_test, rty_test,Fiber_Bending_Energy_std_test,l_p,P0, P1, P2_test, P3, P4, t,boundary_pts,kappa_limit);
                    keep_curve_store(esys_idx)=keep_curve;
                    curve_prob_store(esys_idx)=curve_prob;
                    FiberBendingEnergy_curve_store(esys_idx)=FiberBendingEnergy_curve;
                    
                    if keep_curve>0 || curve_prob==2 % if there is a curve to place or the suggested curve passes outside the boundary
                        % determine which lattice points the proposed curve passes
                        % through and the resulting deformation energy
                        [E_mem,lattice_sat_vals_temp]=Emem_lattice2(rtx,rty,Lat,dA_obstruct_radius,df,lattice_sat_vals);
                        
                        % if the proposed curve passes outside the cell
                        % boundary, compute the resulting deflection energy
                        E_BD=E_BD_curve(rtx_test,rty_test,boundary_pts,t,dist_pair,P0,P1,P2_test,P3,P4,dA_FA_n);
                        
                        % internal energy of the curve, i.e., bending energy
                        Eb=FiberBendingEnergy_curve;
                        
                        % compute proposed energy addition to the system
                        Eprop= Eb+E_mem+E_BD-costs(pt1,pt2)/(2*size(combo_order,1)/numel(unique(combo_order(:))));
                        
                        % compute the energy of the new proposed system
                        E_prop_sys(esys_idx)=E_current_sys_temp(esys_idx)+Eprop;
                        
                        if E_current_sys_temp(esys_idx)>E_prop_sys(esys_idx)
                            place_curve_store(esys_idx)=1;
                            E_new_sys(esys_idx)=E_prop_sys(esys_idx);
                            lattice_sat_vals_test{esys_idx}=lattice_sat_vals_temp;
                        else
                            place_curve_store(esys_idx)=0;
                            E_new_sys(esys_idx)=E_current_sys_temp(esys_idx);
                            lattice_sat_vals_test{esys_idx}=lattice_sat_vals; 
                        end
                    end 
                end % end curve placement determination
                
                % of the 3 potential curves, identify which curves are
                % viable (if any) and record them
                viable_idx1=find((((keep_curve_store>0) + (curve_prob_store==2))>0)==1); % indices of curves that are viable
                viable_idx2=find((place_curve_store>0)==1);
                viable_idx=intersect(viable_idx1,viable_idx2);
                
                if ~isempty(viable_idx)
                    % determine which proposed curve minimizes system
                    % energy. If multiple curves minimize the energy, pick
                    % one at random
                    b_idx=intersect(find(E_new_sys==min(E_new_sys(viable_idx))),viable_idx);
                    if numel(b_idx)>1
                        b_idx2=b_idx(randi(numel(b_idx)));
                    else
                        b_idx2=b_idx;
                    end
                    
                    keep_curve=keep_curve_store(b_idx2);
                    curve_prob=curve_prob_store(b_idx2);
                    place_curve=place_curve_store(b_idx2);
                    rtx=Rtx(b_idx2,:);
                    rty=Rty(b_idx2,:);
                    P2=P2_temp(b_idx2,:);
                    FiberBendingEnergy_curve=FiberBendingEnergy_curve_store(b_idx2);
                    E_new_sys=E_new_sys(b_idx2);
                    lattice_sat_vals=lattice_sat_vals_test{b_idx2};
                else
                    keep_curve=0;
                    curve_prob=0;
                    E_current_sys=E_current_sys_temp(1);
                end

                % RECORD INFORMATION FOR THE CURVE (IF APPLICABLE)
                if keep_curve>0 || curve_prob==2 % if there's a curve to place...
                    if place_curve>0 % if the curve is being placed...
                        % record the curve
                        rt(1,:,bundle_counter)=rtx;
                        rt(2,:,bundle_counter)=rty;

                        % record control points            
                        controlPoints{combo_number}{vec_ctr,bundle_counter}(1,:)=P0;
                        controlPoints{combo_number}{vec_ctr,bundle_counter}(2,:)=P1;
                        controlPoints{combo_number}{vec_ctr,bundle_counter}(3,:)=P2;
                        controlPoints{combo_number}{vec_ctr,bundle_counter}(4,:)=P3;
                        controlPoints{combo_number}{vec_ctr,bundle_counter}(5,:)=P4;

                        % record bending energy                    
                        E_bundle{combo_number}(bundle_counter) = fibers_per_bund*FiberBendingEnergy_curve;

                        % update number of available spots
                        Currentopenpts(pt1)=0; % a 1 indiciates that that lattice point is open; a 0 will indicate that it is occupied
                        Currentopenpts(pt2)=0;

                        % update vector use cell storage
                        W1_temp2=V2_decomp{pt1}; % array of decomposition vectors for point 1
                        W2_temp2=V2_decomp{pt2}; % array of decomposition vectors for point 2
                        I_row_temp=find(ismember(W1_temp2,w1i,'rows')==1);
                        I_col_temp=find(ismember(W2_temp2,w2j,'rows')==1);
                        vecs_used{pt1}(1,I_row_temp) = 1;
                        vecs_used{pt2}(1,I_col_temp) =1;

                        % record time the fibers were created (note: this will
                        % be overwritten as the for loop progresses but it will
                        % be overwritten with the same value so I'm leaving it
                        % here)
                        T_info1=[starting_pt ending_pt I_row_temp I_col_temp];
                        T_info_old=Time_matrix_t{max(tsec_idx-1,1)};
                        if isempty(T_info_old)
                            T_info_old_1=-1.*ones(size(T_info1));
                        else
                            T_info_old_1=T_info_old(:,1:4);
                        end
                        % check if new info exists in old info
                        T_info_check=find(ismember(T_info_old_1,T_info1,'rows')==1);
                        if isempty(T_info_check) % info is new
                            Time_info_new=[T_info1, FiberBendingEnergy_curve, time];
                        else % info is old
                            if bundle_counter<=size(T_info_check,1)
                                Time_info_new=T_info_old(T_info_check(bundle_counter),:);
                                Time_info_new(5)=FiberBendingEnergy_curve;
                            else
                                Time_info_new=[T_info1, FiberBendingEnergy_curve, time];
                            end
                        end

                        Eb_temp_store(vec_ctr,bundle_counter)=FiberBendingEnergy_curve;
                        Lc_temp_store(vec_ctr,bundle_counter)=arclength(rtx,rty) ;
                    else
                        keep_curve=0;
                        Time_info_new=[];
                    end
                else
                    E_new_sys= E_current_sys;
                    Time_info_new=[];
                end

                % record total energy of the system for this pt-pt combo
                E_sys_tot(combo_number)=E_new_sys;

                % increase bundle counter
                bundle_counter=bundle_counter+1;

                % record network time info
                Time_matrix=[Time_matrix; Time_info_new];
            end % end while loop for given P0-P4 combo and the specific vector combo

            % remove rt slices which are all 0, i.e., slices which never describe a curve
            mask = any(any(rt));
            rt = rt(:,:,mask);
            % record rt info for this current vector combo
            rt_vec_storage{vec_ctr}=rt;

            % update W2_temp by removing these vectors from
            % consideration for 
            if ~isempty(W2_temp)
                W2_temp(I_col,:)=[];
            end
        end % end for loop for vector pair determinations

        % remove empty slots from rt_vec_storage
        rt_vec_storage = rt_vec_storage(~cellfun('isempty',rt_vec_storage)); 

        % here's the long way of storing the rt_vec_storage cell into a single
        % 3D matrix:
        % determine total number of fibers created
        fibs=zeros(1,numel(rt_vec_storage)); % vector that'll hold number of fibers made for each vector combo
        for fdet=1:numel(rt_vec_storage)
            ftemp=rt_vec_storage{fdet};
            fibs(fdet)=size(ftemp,3);
        end
        n_fibers_made=sum(fibs); % number of fibers made
        % now just store each of these fibers as a 3D matrix, called rt2
        rt2=zeros(2,length(t),n_fibers_made);
        meh=1;
        for fdet=1:numel(rt_vec_storage)
            ftemp=rt_vec_storage{fdet};
            for fdet2=1:size(ftemp,3)
                rt2(:,:,meh)=ftemp(:,:,fdet2);
                meh=meh+1;
            end
        end

        % record curves for this pt-pt combo
        network{combo_number}=rt2;

        % record system energy for this pt-pt combo
        E_sys_tot(combo_number+1)=E_sys_tot(combo_number);

        % record bending energy for each fiber created for this pt-pt combo.
        % here's one way of doing this: convert the array from a matrix into
        % one long row vector. Then only keep values which are >-1
        Eb_store=Eb_temp_store';
        Eb_store=Eb_store(:)';
        Eb_store=Eb_store(Eb_store>-1);
        Eb_network{starting_pt,ending_pt}=Eb_store;
        Eb_network{ending_pt,starting_pt}=Eb_store;

        Lc_store=Lc_temp_store';
        Lc_store=Lc_store(:)';
        Lc_store=Lc_store(Lc_store>-1);    
        Lc_network{starting_pt,ending_pt}=Lc_store;
        Lc_network{ending_pt,starting_pt}=Lc_store;
    end  % end for loop for rho network constructions

    % Display number of unused vectors
    unused_vec_nums=zeros(1,numberOfPts);
    for vecid=1:numberOfPts
        v3=vecs_used{vecid};
        num_used=nnz(v3);
        num_unused=numel(v3)-num_used;
        unused_vec_nums(vecid)=num_unused;
    end
    %disp('number of unused vectors for each lattice point:')
    %unused_vec_nums
    if sum(unused_vec_nums)<0%sum(unused_vec_nums)>0 % if there are unused vectors 
        % Determine which points were problematic
        pt_remove_vec=network_problem_pts(vecs_used,Lat2,Lat,vec_track,numberOfPts);

        if isempty(pt_remove_vec)
            pt_remove=[];
        else
            % all points in pt_remove_vec are problematic points. first,
            % reduce the number of decomposition vectors for these points
            vec_decomp_numbs(ismember(lat_sat_idx,pt_remove_vec))=vec_decomp_numbs(ismember(lat_sat_idx,pt_remove_vec))-2;
            % now, some points may have been removed because the number of
            % vectors in the decomposition is now <0. 
            removal_idx=find(vec_decomp_numbs(ismember(lat_sat_idx,pt_remove_vec))<0);
            keep_idx=find(vec_decomp_numbs(ismember(lat_sat_idx,pt_remove_vec))>0);
            
            removal_pts=pt_remove_vec(removal_idx);
            altered_points=pt_remove_vec(keep_idx);
            altered_points_vals=vec_decomp_numbs(ismember(lat_sat_idx,altered_points));
            
            % update array of altered points:
            for al=1:numel(altered_points)
                al1=find(NF_altered(:,1)==altered_points(al));
                if isempty(al1)
                    NF_altered=[NF_altered; altered_points(al), altered_points_vals(al) ];
                else
                    NF_altered(al1,2)=altered_points_vals(al);
                end
            end
            
            % if there is a point that needs to be removed, update the
            % pt_remove list
            if ~isempty(removal_idx)
                pt_remove=removal_pts;
            end
            
            pt_remove_store=unique([pt_remove_store,pt_remove]);
        end
    else
        pt_remove=[];
    end
    ncnts=ncnts+1; % loop counter (to avoid infinite loops)
end % end while loop

% remove empty slots from network cell
network2 = network(~cellfun('isempty',network));
    
% record energy of whole system
if size(combo_order,1)>1
    E_sys_tot_plot=[0; E_sys_tot];
else
    E_sys_tot_plot=[0; E_sys_tot'];
end
E_sys_tot_plot(end)=[];

if ~isempty(E_sys_tot_plot)
    E_sys_final=E_sys_tot_plot(end);
else
    E_sys_final=0;
end

end % end function