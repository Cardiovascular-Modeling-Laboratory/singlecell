% extended force equation; test case using f1*L_avg with vhat; vhat is the
% normalized sum of force vectors used to create the fibers connecting any
% 2 points
function [F, F_adh, F_cyto, F_p, F_n] = Force_units_Lavg_linked_v3(R_bound,Npts_t,f_rho,drx,dry,dA,rho_sat,Lc_network,Time_matrix,lat_sat_idx,V2_decomp,Pp_network,Pn_network,fp_tilde,fn_tilde,dA_FA_p,dA_FA_n,rho_0, outside_segs, inside_segs,bdry_mat,drx_norm,dry_norm,Lambda_e)

dA_a=9e-14;
rho_t= dA_a*rho_0/(dA_FA_p-dA_a);
R_t=rho_t/(rho_0+rho_t);

H=heaviside(R_bound-R_t);
H_mat= repmat(H,1,Npts_t)';
H_bound=heaviside(R_bound-R_t);
Inside_segs=sign(bdry_mat+inside_segs);
%% adhesion force
R_array = repmat(R_bound,1,Npts_t)';
sum_fact_x = R_array.*drx; % R_p(r')*[r-r'] in x - will need to sum up each row
sum_fact_y = R_array.*dry; % R_p(r')*[r-r'] in y - will need to sum up each row
F_cyto_x_temp=f_rho.*dA.*R_bound.*sum(sum_fact_x,2);
F_cyto_y_temp=f_rho.*dA.*R_bound.*sum(sum_fact_y,2);

F_adh_sign=sign([F_cyto_x_temp F_cyto_y_temp]);

F_adh=(1.5e+3).*dA_FA_p.*[R_bound R_bound].*F_adh_sign;

F_adh_sign(sum(abs(F_adh_sign),2)==0,:)=1; % if sum(abs(F_adh_sign),2)==0 then R_bound at these points is 0, i.e. there are no bound integrins. so F_adh=0

F_adh_x=F_adh(:,1);
F_adh_y=F_adh(:,2);

%% cytoskeleton force
R_array = repmat(R_bound,1,Npts_t)';
sum_fact_x = R_array.*drx.*H_mat.*Inside_segs; % R_p(r')*[r-r'] in x - will need to sum up each row
sum_fact_y = R_array.*dry.*H_mat.*Inside_segs; % R_p(r')*[r-r'] in y - will need to sum up each row

%fp_t=10^9; % parameter designed to make units work in the concave cell portion
fp_t=2.5486e+7;

F_cyto_x=f_rho.*dA.*R_bound.*H_bound.*sum(sum_fact_x,2)+fp_t.*R_bound.*H_bound.*sum(R_array.*Lambda_e.*drx_norm.*H_mat.*outside_segs.*dA,2);
F_cyto_y=f_rho.*dA.*R_bound.*H_bound.*sum(sum_fact_y,2)+fp_t.*R_bound.*H_bound.*sum(R_array.*Lambda_e.*dry_norm.*H_mat.*outside_segs.*dA,2); %

%% myofibril force terms
Time_matrix_temp=Time_matrix;

if isempty(Time_matrix)  % there's no network to construct
    Pts_actin=[];
else
    Pts_actin=Time_matrix_temp(:,[1,2]);
end

Pts_full=unique(Pts_actin,'rows'); % list of all fiber starting/ending point combinations in the network
Lc_full=cell(Npts_t,Npts_t); % cell array for storing fiber lengths of the full network
Lc_avg=zeros(Npts_t,Npts_t); % storage for average fiber length of any 2 connected points
Nf=zeros(Npts_t,Npts_t); % storage for average fiber length of any 2 connected points

% storage for probability functions
Pp_full=num2cell(zeros(Npts_t,Npts_t));
Pn_full=num2cell(zeros(Npts_t,Npts_t));
Pi_p_full=zeros(Npts_t,Npts_t);
Pi_n_full=zeros(Npts_t,Npts_t);

v_x=zeros(Npts_t,Npts_t); % x component vector of vhat
v_y=zeros(Npts_t,Npts_t); % y component vector of vhat
Ntot=0;
for combo_number=1:size(Pts_full,1)
    starting_pt=Pts_full(combo_number,1);
    ending_pt=Pts_full(combo_number,2);
    
    % construct Lc_avg
    Lc_actin_SE=Lc_network{starting_pt,ending_pt};
    Lc_actin_ES=Lc_network{ending_pt,starting_pt};    
    
    Lc_full{starting_pt,ending_pt}=Lc_actin_SE.*10^(-6); % length of fibers in units of meter
    Lc_full{ending_pt,starting_pt}=Lc_actin_ES.*10^(-6); % length of fibers in units of meter
    
    Lc_avg(starting_pt,ending_pt)= mean(Lc_full{starting_pt,ending_pt});
    Lc_avg(ending_pt,starting_pt)= mean(Lc_full{ending_pt,starting_pt});
    
    % construct Nf
    Nf(starting_pt,ending_pt)=numel(Lc_full{starting_pt,ending_pt});
    Nf(ending_pt,starting_pt)=numel(Lc_full{ending_pt,starting_pt});
    Ntot=Ntot+numel(Lc_full{ending_pt,starting_pt});
    % construct Pi_p and Pi_n
    Pp_full{starting_pt,ending_pt}=Pp_network{starting_pt,ending_pt};
    Pn_full{starting_pt,ending_pt}=Pn_network{starting_pt,ending_pt};
    Pp_full{ending_pt,starting_pt}=Pp_network{ending_pt,starting_pt};
    Pn_full{ending_pt,starting_pt}=Pn_network{ending_pt,starting_pt};
    
    Pi_p_full(starting_pt,ending_pt)=1-prod(1-Pp_full{starting_pt,ending_pt});
    Pi_n_full(starting_pt,ending_pt)=1-prod(1-Pn_full{starting_pt,ending_pt});
    Pi_p_full(ending_pt,starting_pt)=1-prod(1-Pp_full{ending_pt,starting_pt});
    Pi_n_full(ending_pt,starting_pt)=1-prod(1-Pn_full{ending_pt,starting_pt});
    
    % construct vhat_x and vhat_y
    % actin network related vectors:
    T2_temp1a=find(ismember(Time_matrix_temp(:,[1,2]),[starting_pt, ending_pt],'rows')==1);
    if ~isempty(T2_temp1a) 
        % grab vector information for the starting point
        v_idx_SP=unique(Time_matrix_temp(T2_temp1a,3)); % row indices indicating which vectors in the vector decomposition were used to create the fibers of the given combo
        i2_SP=find(starting_pt==lat_sat_idx);
        W_SP = V2_decomp{i2_SP}; % grab the appropriate vector decomposition matrix from the list of decomposition matrices
        W_vecs_actin_SP_1a=W_SP(v_idx_SP,:); % matrix of starting point vectors used to create the fibers

        % grab vector information for the ending point
        v_idx_EP=unique(Time_matrix_temp(T2_temp1a,4)); % row indices indicating which vectors in the vector decomposition were used to create the fibers of the given combo
        i2_EP=find(ending_pt==lat_sat_idx);
        W_EP = V2_decomp{i2_EP}; % grab the appropriate vector decomposition matrix from the list of decomposition matrices
        W_vecs_actin_EP_1a=W_EP(v_idx_EP,:); % matrix of starting point vectors used to create the fibers
    else
        W_vecs_actin_SP_1a=[];
        W_vecs_actin_EP_1a=[];
    end
    
    T2_temp1b=find(ismember(Time_matrix_temp(:,[1,2]),[ending_pt, starting_pt],'rows')==1); % a column vector indicating where the pt-pt combo shows up for actin info
    if ~isempty(T2_temp1b)
        % grab vector information for the starting point
        v_idx_SP=unique(Time_matrix_temp(T2_temp1b,3)); % row indices indicating which vectors in the vector decomposition were used to create the fibers of the given combo
        i2_SP=find(starting_pt==lat_sat_idx);
        W_SP = V2_decomp{i2_SP}; % grab the appropriate vector decomposition matrix from the list of decomposition matrices
        W_vecs_actin_SP_1b=W_SP(v_idx_SP,:); % matrix of starting point vectors used to create the fibers

        % grab vector information for the ending point
        v_idx_EP=unique(Time_matrix_temp(T2_temp1b,4)); % row indices indicating which vectors in the vector decomposition were used to create the fibers of the given combo
        i2_EP=find(ending_pt==lat_sat_idx);
        W_EP = V2_decomp{i2_EP}; % grab the appropriate vector decomposition matrix from the list of decomposition matrices
        W_vecs_actin_EP_1b=W_EP(v_idx_EP,:); % matrix of starting point vectors used to create the fibers
    else
        W_vecs_actin_SP_1b=[];
        W_vecs_actin_EP_1b=[];
    end
    W_vecs_actin_SP=[W_vecs_actin_SP_1a; W_vecs_actin_SP_1b];
    W_vecs_actin_EP=[W_vecs_actin_EP_1a; W_vecs_actin_EP_1b];   
    
    % combine vector information
    W_vecs_SP=W_vecs_actin_SP;
    W_vecs_EP=W_vecs_actin_EP;
    
    % compute average force vector for starting point
    W_vecs_sum_SP=sum(W_vecs_SP,1); % vector sum
    w1_SP=W_vecs_sum_SP./norm(W_vecs_sum_SP); % summation vector
    if ~isempty(w1_SP)
        v_x(starting_pt,ending_pt)=w1_SP(1);
        v_y(starting_pt,ending_pt)=w1_SP(2);
    end 

    % compute average force vector for ending point
    W_vecs_sum_EP=sum(W_vecs_EP,1); % vector sum
    w1_EP=W_vecs_sum_EP./norm(W_vecs_sum_EP); % summation vector
    if ~isempty(w1_EP)
        v_x(ending_pt,starting_pt)=w1_EP(1);
        v_y(ending_pt,starting_pt)=w1_EP(2);
    end 
end

% premyofibril force
%Nf_rand=(0.024-0.17)*rand(size(Nf))+0.17;
if Ntot==0
    Ntot=1;
end
sum_fact_p_x_actin=Nf.*Pi_p_full.*Lc_avg.*drx_norm.*(Nf./Ntot);
sum_fact_p_y_actin=Nf.*Pi_p_full.*Lc_avg.*dry_norm.*(Nf./Ntot);
F_p_x=fp_tilde.*sum(sum_fact_p_x_actin,2).*dA;
F_p_y=fp_tilde.*sum(sum_fact_p_y_actin,2).*dA;

% nascent myofibril force
sum_fact_n_x_actin=Nf.*Pi_n_full.*Lc_avg.*drx_norm.*(Nf./Ntot);
sum_fact_n_y_actin=Nf.*Pi_n_full.*Lc_avg.*dry_norm.*(Nf./Ntot);
F_n_x=fn_tilde.*sum(sum_fact_n_x_actin,2).*dA;
F_n_y=fn_tilde.*sum(sum_fact_n_y_actin,2).*dA;

%% now put together the force equation
F_x=(1./sqrt(sum(abs(F_adh_sign),2))).*F_adh_x+F_cyto_x+F_p_x+F_n_x;
F_y=(1./sqrt(sum(abs(F_adh_sign),2))).*F_adh_y+F_cyto_y+F_p_y+F_n_y;

% total force
F = [F_x F_y];

% record forces
F_cyto=[F_cyto_x F_cyto_y];
F_p=[F_p_x F_p_y];
F_n=[F_n_x F_n_y];

end

