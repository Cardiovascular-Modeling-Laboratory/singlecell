% File name: geometry_change_force_validation_v1.m
%
% Purpose:
% This function calculates the expected force during uniaxial stretching based on
% the stress-strain relationship and compares it to the adhesion force obtained
% from the simulation.
%
% Inputs:
% - E: Young's modulus of the cell (Pa)
% - nu: Poisson's ratio of the cell
% - stretch_factor: The stretch factor applied to the cell (dimensionless)
% - initial_area: The initial area of the cell before stretching (m^2)
% - L_0: The initial length of the cell in the stretch direction (m)
% - F_store: A 3D array storing the force vectors at each time step and point
%            (from your simulation)
% - time_index: The current time index in the simulation
% - dA: The area element associated with each lattice point (m^2)
% - time: The current simulation time (s)
%
% Outputs:
% - F_t: Expected force from the stress-strain relationship (N)
% - F_adhesion: Total adhesion force from the simulation (N)
% - force_difference: Difference between adhesion force and expected force (N)
% - force_ratio: Ratio of adhesion force to expected force (dimensionless)
%
% Usage:
% [F_t, F_adhesion, force_difference, force_ratio] = ...
%     geometry_change_force_validation_v1(E, nu, stretch_factor, initial_area, L_0, ...
%     F_store, time_index, dA, time);

function [F_t, F_adhesion, force_difference, force_ratio] = ...
    geometry_change_force_validation_v1(E, nu, stretch_factor, initial_area, L_0, ...
    F_store, time_index, dA, time)

% Calculate axial strain
epsilon_axial = stretch_factor - 1;

% Calculate transverse strain using Poisson's ratio
epsilon_transverse = -nu * epsilon_axial;

% Original cross-sectional area (perpendicular to stretch direction)
% Assuming initial_area is the area before stretching
% For uniaxial stretch along x-axis, the transverse area is in y-direction
% Since we're in 2D, consider unit thickness in z-direction
A_0 = initial_area / L_0; % Initial width (m)

% Update cross-sectional area
A_t = A_0 * (1 + epsilon_transverse); % Updated width (m)

% Calculate expected force from stress-strain relationship
F_t = A_t * E * epsilon_axial; % Force in Newtons (N)

% Extract adhesion force from simulation
% Assuming stretch is along x-axis
F_x = squeeze(F_store(time_index, :, 1)); % x-component of force at current time
% Total adhesion force in x-direction
F_adhesion = sum(F_x .* dA); % Sum over all points (N)

% Compare forces
force_difference = F_adhesion - F_t;
force_ratio = F_adhesion / F_t;

% Output results
disp(['Time: ', num2str(time)]);
disp(['Expected Force F_t: ', num2str(F_t), ' N']);
disp(['Adhesion Force F_adhesion: ', num2str(F_adhesion), ' N']);
disp(['Force Difference: ', num2str(force_difference), ' N']);
disp(['Force Ratio: ', num2str(force_ratio)]);

end
