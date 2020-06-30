% File name: R_factors.m
% 
% Last updated: 9/11/2018
% First version finished: 08/21/2009
% 
% written by Anya Grosberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Harvard University, Cambridge, MA 02138
% 
% The purpose of this function is to calculate the R factors given the
% integrin concentrations
% 
% 
% Input:  1. rho_0 -- saturation integrin constant (obtained from
%         parameter_gen_units file
%         2. Pre-myofibril bound integrin concentrations
%         3. Nascent-myofibril bound integrin cocentration
%        
%         
% Output: 1. R_n (Look at equations for definition)
%         2. R_p (Look at equations for definition)


function R_bound = R_factors_units_v2(rho_0,integrin_bound)

Sat_denom = rho_0+integrin_bound; %the denomenator, which is the same for both R factors
R_bound = integrin_bound./Sat_denom; %R factor for the pre-myofibril integrin
end
