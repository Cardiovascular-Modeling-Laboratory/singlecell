% free energy model for second order phase transitions modeled after
% soukhovolsky et al. A Boltzmann-type probability is used to estimate the
% probability of a curve being in state 1 (premyofibril) or state 2
% (nascent myofibril) based on an energy function. The energy is the
% "maturation" energy where the energy of state 1 equals the bending energy
% and the energy of state 2 follows the energy described in soukhovolsky et
% al.

function [Pp_network, Pn_network,fiber_id1] = fiber_id_linked2(Time_matrix_t,Npts_t, Eb_max, fibers_per_bund,time_index,a,b,F_x, F_y, Fsat_lat, F_c)

ti=time_index+1;
Time_matrix=Time_matrix_t{ti};

Pp_network=num2cell(zeros(Npts_t,Npts_t));
Pn_network=num2cell(zeros(Npts_t,Npts_t));
fiber_id1=cell(Npts_t,Npts_t);

% magnitude of force at each lattice point:
Fmag=sqrt(F_x.^2+F_y.^2);

% determine probability function for the actin network
for T_el=1:size(Time_matrix,1)
    starting_pt=Time_matrix(T_el,1);
    ending_pt=Time_matrix(T_el,2);
    EB=Time_matrix(T_el,5);
%    Tbuild=Time_matrix(T_el,6);
    F_sp=Fmag(starting_pt);
    F_ep=Fmag(ending_pt);
    F_min = min(F_sp, F_ep); % minimal force being generated at each endpoint
   
    E1=EB;
    if F_min>=Fsat_lat && F_min<=F_c
        E2=Eb_max+3;
    else
        E2=Eb_max+3-((a^2)/(4*b))*((1./(F_min-Fsat_lat)-1./F_c).^2);
    end

    deltaE=E2-E1;
    P_p=1./(1+exp(-deltaE));
    P_n=1./(1+exp(deltaE));
    
    Pp_network{starting_pt,ending_pt}=[Pp_network{starting_pt,ending_pt} P_p];
    Pn_network{starting_pt,ending_pt}=[Pn_network{starting_pt,ending_pt} P_n];
    Pp_network{ending_pt,starting_pt}=[Pp_network{ending_pt,starting_pt} P_p];
    Pn_network{ending_pt,starting_pt}=[Pn_network{ending_pt,starting_pt} P_n];

    % determine the state of each costructed fiber
    prob=[P_p P_n];
    if sum(prob)>0 % if at least one probability is non-zero
        samp=randsample(1:numel(prob), fibers_per_bund, true, prob);
        n1=sum(samp==1); % determine how many times "1" was chosen (premyofibril designation)
        n2=sum(samp==2); % determine how many times "2" was chosen (nascent myofibril designation)
        [~,max_idx]=max([n1, n2]);
        id_temp=max_idx;
    else
        id_temp=[];
    end
    
    fiber_id1{starting_pt,ending_pt}=[fiber_id1{starting_pt,ending_pt} id_temp];
    fiber_id1{ending_pt,starting_pt}=[fiber_id1{ending_pt,starting_pt} id_temp];

end

end % end function



