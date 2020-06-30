% File name: parameter_gen_units_linked_v2.m
%
%
% Last updated: 4/15/2020 (by Will)
% First version finished: 08/31/2009
% 
% written by Anya Grosberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Harvard University, Cambridge, MA 02138
% 
% This code creates a .mat file that contains all of the parameter values
% used throughout the simulation. The user can edit the parameter values as
% they wish, although they really shouldn't edit many of the parameter
% values since they were obtained through fitting and don't change when
% cell size/shape are changed.
% The user can save the parameter file (the last line of the code) however
% they want. I advise that it's saved in a consistent manner or using some
% kind of summary file to keep track of what differentiates each saved
% parameter file.
% For reference/more info, look at the Tables for model implementation in
% my notes/write up
% 
% Inputs:  None; simply Run this code before proceeding with the simulation
%        
%         
% Output: .mat file with all the constants that will be used in the
% simulation

function parameter_gen_units_linked_v2()
%% Fixed constants: These constants SHOULD NOT be changed
F_0=3.8427e-09; % ~4nN in units of newtons

% ----- ----- ----- ----- Force equation related
Fsat_lat=3.0e-9; % force limit required to begin construction premyofibrils
F_sat_n=9e-9; % force limit required for premyofibril to mature into nascent myofibril
dA_FA_n=8*10^(-12); % area of focal adhesion for nascent myofibril (converted from um^2 to m^2)
dA_FA_p=2*10^(-12); % area of focal adhesion for premyofibril (converted from um^2 to m^2)

% ----- ----- ----- ----- R related
rho_0 = 2.8e+14; %3.08e+13; %the saturation value of the bound integrin cocnetration
R_sat = 0.9; % parameter used to define rho_sat 0.91

% ----- ----- ----- ----- Fiber model related (in micrometer scale since
% the fiber model codes are written in the micrometer scale!!!)
rc_limit=0.3; % limit of radius of curvature at any point on a constructed curve (in um). must be >=0.3
kappa_limit=1/rc_limit;
df=0.14; % diameter of myofibril bundle in micometer: can be between (0.133,0.306) um
fibers_per_bund=20;
l_p = 17; % persistence length of actin in micrometers
t0=1.5;%1.3841; % thickness of a bundle in micrometers
Eb_max= 13.1; % max "noramlized" bending energy, in units of 1/um: see code notes for where this came from
fib_sat_limit=1; % record fiber saturation limit the lattice %Rsat_lat=0.377; % this is the saturation limit used to determine if a lattice point is sufficiently saturated

max_bundle_number = 3; %maximum number of curves to attempt to create during fiber creation 
kt=15;

% ----- ----- ----- ----- Fiber model type identification related

a=2*51.1; % second order phase transition parameter (in J*s)
Fmax=18e-9; % max force after which the fiber has almost certainly matured
Fm = Fmax-Fsat_lat; % the additional force needed to get from F^* (Fsat_lat in the code) to Fmax: Fmax = F_m+F^*
F_c=  F_sat_n; % critical force threshold needed to begin maturation 
b=(a^2/4)*(1/(Eb_max+3+log(99)))*(1/Fm-1/F_c)^2;

% ----- ----- ----- ----- Simulation related (general)
rho_bar = 400e+12; % From literature reference (see notes)
t_hr=72; % desired simulation time (in hrs) 

t_max_u=3600*t_hr; % implemented simulation time (in sec) 
t_max=(120/259200)*t_max_u; % dimensionless simulation time (in arbitrary units)
time_step = 0.08; % delta t
time_max_factor = floor(t_max/time_step); % this defines the number of non-zero time points, N_t in the notes
time_step_u = t_max_u/time_max_factor; % delta t^u

% ----- ----- ----- ----- Integrin rate constants 
k_0=1.6223e-05;
k_1=3.2640e+03;
k_m1=2.8063e-04;

%% Flexible constants: These constants CAN be changed from simulation to simulation
A = 2116*10^(-12); % Area of the cell written here in square meters, |Delta| in the notes
Tot_Integrin_Val = rho_bar*A; % Total amount of integrin with units included, N_p in the notes
rho_sat = R_sat*rho_0*(1/(1-R_sat)); % max denisty of integrins per unit area of the cell

Num_points=100; % The minimum number of points in the grid, N_pts in the notes

nuc_area = 100*10^(-12); % The area of the nucleus in m for most cells
%nuc_area = 50*10^(-12); % The area of the nucleus in m for most 1250um^2 cell. It's smaller so that I can run the code for elongated cells and nucleus is contained in the elongated cell
nuc_radius=sqrt(nuc_area/pi);

% ----- ----- ----- ----- Color map
% Save the green colormap for the bound integrins, blue colormap for the
% unbound integrin
bound_colormap = [0 0.2 0.2;0 0.2127 0.1968;0 0.2254 0.1937;0 0.2381 0.1905;0 0.2508 0.1873;0 0.2635 0.1841;0 0.2762 0.181;0 0.2889 0.1778;0 0.3016 0.1746;0 0.3143 0.1714;0 0.327 0.1683;0 0.3397 0.1651;0 0.3524 0.1619;0 0.3651 0.1587;0 0.3778 0.1556;0 0.3905 0.1524;0 0.4032 0.1492;0 0.4159 0.146;0 0.4286 0.1429;0 0.4413 0.1397;0 0.454 0.1365;0 0.4667 0.1333;0 0.4794 0.1302;0 0.4921 0.127;0 0.5048 0.1238;0 0.5175 0.1206;0 0.5302 0.1175;0 0.5429 0.1143;0 0.5556 0.1111;0 0.5683 0.1079;0 0.581 0.1048;0 0.5937 0.1016;0 0.6063 0.09841;0 0.619 0.09524;0 0.6317 0.09206;0 0.6444 0.08889;0 0.6571 0.08571;0 0.6698 0.08254;0 0.6825 0.07937;0 0.6952 0.07619;0 0.7079 0.07302;0 0.7206 0.06984;0 0.7333 0.06667;0 0.746 0.06349;0 0.7587 0.06032;0 0.7714 0.05714;0 0.7841 0.05397;0 0.7968 0.05079;0 0.8095 0.04762;0 0.8222 0.04444;0 0.8349 0.04127;0 0.8476 0.0381;0 0.8603 0.03492;0 0.873 0.03175;0 0.8857 0.02857;0 0.8984 0.0254;0 0.9111 0.02222;0 0.9238 0.01905;0 0.9365 0.01587;0 0.9492 0.0127;0 0.9619 0.009524;0 0.9746 0.006349;0 0.9873 0.003175;0 1 0];
unbound_colormap = [0 0 0.2;0 0.01587 0.2127;0 0.03175 0.2254;0 0.04762 0.2381;0 0.06349 0.2508;0 0.07937 0.2635;0 0.09524 0.2762;0 0.1111 0.2889;0 0.127 0.3016;0 0.1429 0.3143;0 0.1587 0.327;0 0.1746 0.3397;0 0.1905 0.3524;0 0.2063 0.3651;0 0.2222 0.3778;0 0.2381 0.3905;0 0.254 0.4032;0 0.2698 0.4159;0 0.2857 0.4286;0 0.3016 0.4413;0 0.3175 0.454;0 0.3333 0.4667;0 0.3492 0.4794;0 0.3651 0.4921;0 0.381 0.5048;0 0.3968 0.5175;0 0.4127 0.5302;0 0.4286 0.5429;0 0.4444 0.5556;0 0.4603 0.5683;0 0.4762 0.581;0 0.4921 0.5937;0 0.5079 0.6063;0 0.5238 0.619;0 0.5397 0.6317;0 0.5556 0.6444;0 0.5714 0.6571;0 0.5873 0.6698;0 0.6032 0.6825;0 0.619 0.6952;0 0.6349 0.7079;0 0.6508 0.7206;0 0.6667 0.7333;0 0.6825 0.746;0 0.6984 0.7587;0 0.7143 0.7714;0 0.7302 0.7841;0 0.746 0.7968;0 0.7619 0.8095;0 0.7778 0.8222;0 0.7937 0.8349;0 0.8095 0.8476;0 0.8254 0.8603;0 0.8413 0.873;0 0.8571 0.8857;0 0.873 0.8984;0 0.8889 0.9111;0 0.9048 0.9238;0 0.9206 0.9365;0 0.9365 0.9492;0 0.9524 0.9619;0 0.9683 0.9746;0 0.9841 0.9873;0 1 1];
fiber_colormap = [1,1,1;0.976470589637756,0.976470589637756,0.976470589637756;0.952941179275513,0.952941179275513,0.952941179275513;0.929411768913269,0.929411768913269,0.929411768913269;0.905882358551025,0.905882358551025,0.905882358551025;0.870898067951202,0.870979607105255,0.870898067951202;0.835913777351379,0.836076796054840,0.835913777351379;0.800929427146912,0.801174044609070,0.800929427146912;0.765945136547089,0.766271233558655,0.765945136547089;0.730960845947266,0.731368482112885,0.730960845947266;0.695976555347443,0.696465730667114,0.695976555347443;0.660992205142975,0.661562919616699,0.660992205142975;0.626007914543152,0.626660168170929,0.626007914543152;0.591023623943329,0.591757416725159,0.591023623943329;0.556039333343506,0.556854605674744,0.556039333343506;0.521054983139038,0.521951854228973,0.521054983139038;0.486070692539215,0.487049072980881,0.486070692539215;0.451086401939392,0.452146291732788,0.451086401939392;0.443856894969940,0.444860994815826,0.443856894969940;0.436627358198166,0.437575668096542,0.436627358198166;0.429397851228714,0.430290371179581,0.429397851228714;0.422168314456940,0.423005074262619,0.422168314456940;0.414938807487488,0.415719777345657,0.414938807487488;0.407709270715714,0.408434450626373,0.407709270715714;0.400479763746262,0.401149153709412,0.400479763746262;0.393250226974487,0.393863856792450,0.393250226974487;0.386020720005035,0.386578559875488,0.386020720005035;0.378791183233261,0.379293233156204,0.378791183233261;0.371561676263809,0.372007936239243,0.371561676263809;0.364332139492035,0.364722639322281,0.364332139492035;0.357102632522583,0.357437342405319,0.357102632522583;0.349873095750809,0.350152015686035,0.349873095750809;0.342643588781357,0.342866718769074,0.342643588781357;0.335414052009583,0.335581421852112,0.335414052009583;0.328184545040131,0.328296124935150,0.328184545040131;0.320955008268356,0.321010798215866,0.320955008268356;0.313725501298904,0.313725501298904,0.313725501298904;0.309111893177032,0.309111893177032,0.309111893177032;0.304498285055161,0.304498285055161,0.304498285055161;0.299884676933289,0.299884676933289,0.299884676933289;0.295271068811417,0.295271068811417,0.295271068811417;0.290657460689545,0.290657460689545,0.290657460689545;0.286043822765350,0.286043822765350,0.286043822765350;0.281430214643478,0.281430214643478,0.281430214643478;0.276816606521606,0.276816606521606,0.276816606521606;0.272202998399735,0.272202998399735,0.272202998399735;0.267589390277863,0.267589390277863,0.267589390277863;0.262975782155991,0.262975782155991,0.262975782155991;0.258362174034119,0.258362174034119,0.258362174034119;0.253748565912247,0.253748565912247,0.253748565912247;0.249134957790375,0.249134957790375,0.249134957790375;0.244521334767342,0.244521334767342,0.244521334767342;0.239907726645470,0.239907726645470,0.239907726645470;0.235294118523598,0.235294118523598,0.235294118523598;0.211764708161354,0.211764708161354,0.211764708161354;0.188235297799110,0.188235297799110,0.188235297799110;0.164705887436867,0.164705887436867,0.164705887436867;0.141176477074623,0.141176477074623,0.141176477074623;0.117647059261799,0.117647059261799,0.117647059261799;0.0941176488995552,0.0941176488995552,0.0941176488995552;0.0705882385373116,0.0705882385373116,0.0705882385373116;0.0470588244497776,0.0470588244497776,0.0470588244497776;0.0235294122248888,0.0235294122248888,0.0235294122248888;0,0,0];

filename = ['paramTEST5_' num2str(t_hr) 'hr.mat'];

save(filename,'a','A',...
    'b','bound_colormap',...
    'dA_FA_p','dA_FA_n','df',...
    'Eb_max',...
    'F_0','F_sat_n','Fsat_lat','fib_sat_limit','fibers_per_bund','Fmax','Fm','F_c','fiber_colormap',...
    'kappa_limit','kt','k_0','k_1','k_m1',...
    'l_p',...
    'max_bundle_number',...
    'Num_points','nuc_area','nuc_radius',...
    'R_sat','rho_0','rc_limit','rho_bar','rho_sat',...
    't0','t_hr','t_max_u','t_max','time_step','time_max_factor','time_step_u','Tot_Integrin_Val',...
    'unbound_colormap');
