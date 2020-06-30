function [E_prop_sys_nesp_temp,lattice_sat_vals_temp2,rt_temp2,ep1,pt_remove_vec2,v2_store_x,v2_store_y,Lc_n,Eb_n]=forced_curve(P0,v1,Lat2,t,lattice_sat_vals_temp,kappa_limit,l_p,boundary_pts,E_current_sys,Lat,...
    dA_obstruct_radius,kt,df,brack,lattice_sat_vals,V2,Eb_max)

pt_remove_vec2=[];
P0_idx=find(ismember(Lat2,P0,'rows')==1);

% possible ending points are all other saturated lattice points
ending_pt_vec=Lat2;
ending_pt_vec(P0_idx,:)=[];
V2_new=V2;
V2_new(P0_idx,:)=[];

E_prop_sys_nesp=10^5.*ones(1,size(ending_pt_vec,1));
rt_nesp=cell(1,size(ending_pt_vec,1));
lattice_sat_vals_temp_matrix=repmat(lattice_sat_vals_temp,1,size(ending_pt_vec,1));
Eb_store=zeros(1,size(ending_pt_vec,1));
Lc_store=zeros(1,size(ending_pt_vec,1));
keep_curve_vec=zeros(2,size(ending_pt_vec,1));
v2_store=zeros(size(ending_pt_vec,1),2);

for end_pt_idx=1:size(ending_pt_vec,1)
    rt_temp=zeros(2,length(t));
    ending_pt=ending_pt_vec(end_pt_idx,:);
    P4=ending_pt;
    
    v2=V2_new(end_pt_idx,:);
%     % v2 and v1 may have norms which are vastly different in
%     % magnitude. To avoid this being a potential issue, rescale v2
%     % so that it is of a similar magnitude to v1. To do this,
%     % determine the relative scale of v1 and then renormalize v2 to
%     % be on the same scale
%     v1_temp=v1;
%     v1_temp(v1_temp==0)=[];
%     cv1=max(floor(log10(abs(v1_temp))));
%     v2=(v2./norm(v2)).*10^cv1; 
     v2_store(end_pt_idx,:)=v2;
    
    % construct min energy curve, assuming curve avoids obstructions
    [P1,P2,P3,Fiber_Bending_Energy_std,rtx,rty] = MinEnergyCurveWithoutObstruct3(P0,P4,v1,v2,t,kappa_limit,0,[],0,0,0);
    
    [FiberBendingEnergy_curve, keep_curve,curve_prob]=curve_comps(rtx, rty,Fiber_Bending_Energy_std,l_p,P0, P1, P2, P3, P4, t,boundary_pts,kappa_limit); 
    keep_curve_vec(1,end_pt_idx)=keep_curve;
    keep_curve_vec(2,end_pt_idx)=curve_prob;
    
    % CURVE PLACEMENT VIABLE?
    if keep_curve>0  % if there is a curve to place within the cell
        % determine which lattice points the proposed curve passes
        % through and the resulting deformation energy
        [E_mem,lattice_sat_vals_temp]=Emem_lattice(rtx,rty,Lat,dA_obstruct_radius,kt,df,brack,lattice_sat_vals);
        lattice_sat_vals_temp_matrix(:,end_pt_idx)=lattice_sat_vals_temp;

        % internal energy of the curve, i.e., bending energy
        Eb=FiberBendingEnergy_curve;
        Eb_store(end_pt_idx)=Eb;

        % compute proposed energy addition to the system
        Eprop= (Eb-Eb_max)+E_mem;

        % compute the energy of the new proposed system
        E_prop_sys_nesp(end_pt_idx)=E_current_sys+Eprop;

        % temporarily record the curve 
        rt_temp(1,:)=rtx;
        rt_temp(2,:)=rty;
        rt_nesp{end_pt_idx}=rt_temp;

        % record length of proposed curve
        Lc_store(end_pt_idx)=arclength(rtx,rty);
    else
        rt_nesp{end_pt_idx}=rt_temp;
    end % end curve placement determination
    
end

% place best curve
curve0=find(keep_curve_vec(1,:)==1); % if there's an index with keep_curve=1, grab the best one
if ~isempty(curve0)
    E_prop_sys_nesp_c0=E_prop_sys_nesp(curve0); 
    [E_prop_sys_nesp_temp,mi2]=min(E_prop_sys_nesp_c0);
    E_prop_sys_nesp_min_idx=curve0(mi2);            
    v2_store_x=v2_store(E_prop_sys_nesp_min_idx,1);
    v2_store_y=v2_store(E_prop_sys_nesp_min_idx,2);
    Lc_n=Lc_store(E_prop_sys_nesp_min_idx);
    Eb_n=Eb_store(E_prop_sys_nesp_min_idx);    
else
    % every keep_curve is 0. So the point needs to be removed from
    % consideration.
    % the next few lines are superfluous. they're just there so
    % that i can proceed with the code without errors but will not play a
    % role in the long run. all i care about at this point is that the
    % lattice point of interest shouldn't be used to create the network. so
    % record it in pt_remove_vec2
    E_prop_sys_nesp_min_idx=1;
    E_prop_sys_nesp_temp=E_prop_sys_nesp(E_prop_sys_nesp_min_idx); 
    
    v2_store_x=v2_store(E_prop_sys_nesp_min_idx,1);
    v2_store_y=v2_store(E_prop_sys_nesp_min_idx,2);
    Lc_n=Lc_store(E_prop_sys_nesp_min_idx);
    Eb_n=Eb_store(E_prop_sys_nesp_min_idx);
    
    P0_idx2=find(ismember(Lat,P0,'rows')==1);
    pt_remove_vec2=[pt_remove_vec2,P0_idx2];
end

lattice_sat_vals_temp2=lattice_sat_vals_temp_matrix(:,E_prop_sys_nesp_min_idx); % update lattice saturation values

rt_temp2=rt_nesp{E_prop_sys_nesp_min_idx};
ep1=ending_pt_vec(E_prop_sys_nesp_min_idx,:);

end % end function