% This code creates a template network and is used to determine the order
% in which the network will be constructed (combo_order)

function [A_full, A_template, combo_order]=network_template(Lat,outline, V2, Lat2,outside_segs, inside_segs,bdry_mat,lat_sat_idx,choice,vec_decomp_numbs)

if isempty(Lat2)
%    disp('  Error: No network to be constructed. Try increasing energy benefit')
%    emailTest();
    A_full=[]; A_template=[]; combo_order=[];
    return
end

% graph information
Npts=size(Lat2,1);
r_col_x = repmat(Lat2(:,1),1,Npts); %make the array with each column having rx coorinates
r_row_x = r_col_x'; %make the array with each row having rx coorinates
r_col_y = repmat(Lat2(:,2),1,Npts); %make the array with each column having ry coorinates
r_row_y = r_col_y';
drx = r_row_x - r_col_x; %the differences in x coordinates - each row - vary r', column - vary r
dry = r_row_y - r_col_y;
dr_dist_squared = drx.^2+dry.^2; %the distance between each r'-r pair squared

L_ij = sqrt(dr_dist_squared);
diag_indices=linspace(1,Npts^2,Npts); 
L_ij(diag_indices)=0;

% construct theta_ij and delta_s. Note that these should be symmetric in
% (i,j) so construct them as triangular matrices
% Since we're requiring a_ij=0 in the simulation, in order for this case to
% not contribute to the objective function, we can set theta_ij=pi/2 and
% delta_s=0
theta_ij_temp=zeros(Npts,Npts);

for sim_idx=1:Npts
    P0=Lat2(sim_idx,:);    
    W1=V2(sim_idx,:);
    for sim_idx2=sim_idx:Npts
        P4=Lat2(sim_idx2,:);        
        W2=V2(sim_idx2,:);
        [theta_ij_val, ~]=vector_similarity2(P0,P4,W1,W2,outline);
        theta_ij_temp(sim_idx,sim_idx2)=theta_ij_val;
%        delta_s_temp(sim_idx,sim_idx2)=delta_ij_val;
    end
end
theta_ij=theta_ij_temp+theta_ij_temp';
theta_ij(diag_indices)=pi/2;
%delta_s=delta_s_temp+delta_s_temp';
delta_s_temp=sign(bdry_mat+inside_segs+outside_segs); % all point-combos that can be considered (using Lat)

L_ij_temp=zeros(size(L_ij));
L_ij_temp(L_ij>25)=1;

delta_s=delta_s_temp(lat_sat_idx,lat_sat_idx).*L_ij_temp;% all point-combos that can be considered (using Lat2)

%delta_s=delta_s_temp(lat_sat_idx,lat_sat_idx);
%delta_s=ones(Npts,Npts);
delta_s(diag_indices)=0;


% create complete graph
A_full=ones(Npts,Npts);
A_full(diag_indices)=0;
A_full=A_full.*delta_s; 
%% determine filtered graph
% so far, we have theta_ij which determines how similar the 2 vectors are
% if they were placed at the same starting point. but it doesn't take into
% account the starting positions of each vector, i.e., whether 2 vectors
% point "towards each other" or point "away from each other". we need to
% determine this in order to build the appropriate template network.

%c_ij_temp=sign(bdry_mat+inside_segs)-sign(outside_segs);
%c_ij=abs(c_ij_temp(lat_sat_idx,lat_sat_idx)); % matrix of pm1 depending on whether the connection is convex (1), concave (-1), or n/a (0)
c_ij_temp=sign(bdry_mat+inside_segs)+sign(outside_segs);
c_ij=c_ij_temp(lat_sat_idx,lat_sat_idx); % matrix of pm1 depending on whether the connection is convex (1), concave (-1), or n/a (0)

% determine the points which are along the concave region of the cell
concave_idx=unique(find(sum(outside_segs(lat_sat_idx,lat_sat_idx))>0));
replace_idx=setdiff(1:size(c_ij,1),concave_idx);
c_ij_concave = -1.*abs(c_ij(concave_idx', replace_idx'));
%c_ij_old=c_ij;
c_ij_new=c_ij;
c_ij_new(concave_idx',replace_idx')=c_ij_concave;
c_ij_new(replace_idx',concave_idx')=c_ij_concave';
c_ij=c_ij_new;

H_A=zeros(1,Npts^2);

A_old=A_full;
H_A(1)=0.5*sum(sum(c_ij.*A_old.*theta_ij));

for i=2:Npts^2
    A_test=A_old;
    [i_sub,j_sub]=ind2sub([Npts,Npts],i);
    A_test(i_sub,j_sub)=0;
    A_test(j_sub,i_sub)=0;
    H_A_old=0.5*sum(sum(c_ij.*A_old.*theta_ij));
    H_A_test=0.5*sum(sum(c_ij.*A_test.*theta_ij));
    
    if H_A_old< H_A_test
        A_new=A_old;
    else
        A_new=A_test;
    end
    
    H_A(i)=0.5*sum(sum(c_ij.*A_new.*theta_ij));
    A_old=A_new;
end

A_filtered=A_old;

Vmag=sqrt(V2(:,1).^2+V2(:,2).^2);
[~,sorting_idx]=sort(Vmag);
sorted_indices=fliplr(sorting_idx')';


A_filtered2=sign(A_filtered+A_filtered');

% determine template network
Atemp=zeros(size(A_filtered2));

%Mfibs2=min(vec_decomp_numbs);
comb_val=[];

for i=1:numel(sorted_indices)
    id1=sorted_indices(i);
    Mfibs2=vec_decomp_numbs(i); %_net2
    m2=setdiff(find(A_filtered2(id1,:)==1),sorted_indices(1:i)); % list of potential connecting points for point id1
    
    if ~isempty(m2)
        [~, ~, m2_sorted] = m2sort(id1, Lat2, V2, m2,theta_ij,outside_segs,dr_dist_squared);
        i_rel=m2_sorted(1:min(min(Mfibs2,Mfibs2-sum(comb_val(:)==id1)),numel(m2_sorted)))';%Mfibs2-(numel(find(A_filtered2(id1,:)==1))-numel(m2)))))';
    else
        i_rel=[];
    end
    
    comb_vals_temp=[id1.*ones(size(i_rel)), i_rel];
    comb_val=[comb_val; comb_vals_temp];
end

for i=1:size(comb_val,1)
    Atemp(comb_val(i,1),comb_val(i,2))=1;
    Atemp(comb_val(i,2),comb_val(i,1))=1;
end
A_template=Atemp;

%% construct order of network construction, combo_order
combo_order_temp=[];
for j=1:Npts
    sp=sorted_indices(j);
    
    for j2=1:Npts
        eps=find(A_template(sp,:)==1);
        % it's possible that what's best for the system is if no connection
        % is made to a particular point. in this case, the point should not
        % be considered in the combo order; a code later in the fiber model
        % will remove points that fall into this category
        if ~isempty(eps)
            [costs,paths] = dijkstra(A_template,L_ij,sp,eps);
            [~,sorted_cost_idx]=sort(costs);
            if numel(sorted_cost_idx)>1
                sorted_paths=paths(sorted_cost_idx);
            else
                sorted_paths=cell(1,1);
                sorted_paths{1}=paths;
            end

            for j3=1:numel(sorted_cost_idx)
                shortest_path=sorted_paths{j3};
                for comb_id=1:max(numel(shortest_path)-1,1)
                    combo_order_temp=[combo_order_temp; shortest_path(comb_id), shortest_path(comb_id+1)];
                end            
            end
        end
    end
    

end

% remove duplicate rows from combo_order (with Lat2 indices)
[~,id2]=ismember(unique(sort(combo_order_temp,2),'rows'),sort(combo_order_temp,2),'rows');
combo_order_Lat2=combo_order_temp(sort(id2),:); % combo order using Lat2 indexing

% convert combo_order from Lat2 indexing to Lat indexing
Lat2_idxing=[];
Lat_idxing=[];
for j=1:Npts
    P0=Lat2(j,:);
    Lat_idxing=[Lat_idxing; j find(ismember(Lat,P0,'rows')==1)];
    Lat2_idxing=[Lat2_idxing; j.*ones(size(find(combo_order_Lat2==j))) find(combo_order_Lat2==j)];
end

replacement=zeros(size(Lat2_idxing));
for j=1:size(replacement,1)
    replacement(j,1)=Lat2_idxing(j,2);
    k=Lat2_idxing(j,1);
    replacement(j,2)=Lat_idxing(find(Lat_idxing(:,1)==k),2);
end
combo_order=combo_order_Lat2;
if ~isempty(combo_order)
    combo_order(replacement(:,1))=replacement(:,2);
    % due to the construction approach, we will always require the fibers to be
    % constructed from lattice point i to lattice point j where i<j. Hence, the
    % first element of each row of combo_order must be smaller than the 2nd
    % element of that row:
    combo_order=sort(combo_order,2);
end

%% plot the template network
% Npts=size(Lat2,1);
% SP_set=[];
% EP_set=[];
% xplot=reshape(A_template,Npts,Npts);
% %xplot=x0;
% 
% for plot_idx1=1:Npts
%     EP=Lat2(find(xplot(plot_idx1,1:Npts)==1),:);
%     
%     SP=repmat(Lat2(plot_idx1,:),size(EP,1),1);
%         
%     SP_set=[SP_set; SP];
%     EP_set=[EP_set; EP];
% end
% 
% figure
% plot(Lat2(:,1),Lat2(:,2),'ko')
% hold on
% plot(Lat(:,1),Lat(:,2),'k.')
% for j=1:size(SP_set,1)
%     line_set=[SP_set(j,:); EP_set(j,:)];
%     line(line_set(:,1),line_set(:,2))
% end

end % end function