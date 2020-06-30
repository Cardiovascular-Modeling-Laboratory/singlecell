% File name: New_integrin_cons_units.m
% 
% Last updated: 9/12/2018
% First version finished: 08/25/2009
% 
% written by Anya Grosberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Harvard University, Cambridge, MA 02138
% 
% The purpose of this function is to calculate the concentrations of the
% three types of integrins at the next time step
% 
% 
% Input:  1. integrin_free - concentration of the free integrin at time t_1
%         2. integrin_p - concentration of the bound integrin connected to pre-myofibril at time t_1
%         3. integrin_n - concentration of the bound integrin connected to nascent myofibrils at time t_1
%         4. F - net force, from Force_units.m
%         5. U - bias potential, from Bias_Pot_units.m
%      6-11. k_0,k_1,k_m1,k_2,k_m2,tau - constants from parameter_gen_units
%      file
%        12. time_step_u - the step of time that we are taking, from
%        parameter_gen_units
%        13. rho_bar - average integrin density in the cell, from
%        parameter_gen_units file
%        14. Concave_ind - a vector the same length as mat_r, but with
%        zeros where there is no ECM in the concave shapes
%        
%         
% Output: 1. integrin_free_new - concentration of the free integrin at time t_2
%         2. integrin_p_new - concentration of the bound integrin connected to pre-myofibril at time t_2
%         3. integrin_n_new  - concentration of the bound integrin connected to nascent myofibrils at time t_2


% function [integrin_free_new, integrin_bound_new]=New_integrin_cons_units_v2(integrin_free, integrin_bound,...
%     F, U,k_0,k_1,k_m1,k_2,tau,time_step_u,rho_bar,Concave_ind,F_0)
function [integrin_free_new, integrin_bound_new]=New_integrin_cons_units_v2(integrin_free, integrin_bound,...
    F, k_0,k_1,k_m1,time_step_u,rho_bar,Concave_ind,F_0)
%global Fig_Num

%find the magnitude of F at each r
absF = sqrt(sum((F.^2),2));

%the myofibril connected integrins
integrin_bound_new = ((k_0+k_1.*absF).*integrin_free -k_m1.*exp(-absF./F_0) .*integrin_bound ).*time_step_u+integrin_bound;

%free integrins
integrin_free_new = Concave_ind.*((sum((rho_bar-integrin_bound_new))./sum(Concave_ind.*ones(size(integrin_bound_new)))));
