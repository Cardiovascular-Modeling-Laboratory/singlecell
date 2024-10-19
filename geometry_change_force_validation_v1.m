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

% Convert units where necessary
initial_area = initial_area * 1e-12; % Convert from µm² to m²
dA = dA * 1e-12;                     % Convert from µm² to m²

% Calculate axial strain
epsilon_axial = stretch_factor - 1;

% Calculate transverse strain using Poisson's ratio
epsilon_transverse = -nu * epsilon_axial;

% Original cross-sectional width (perpendicular to stretch direction)
A_0 = initial_area / L_0; % Initial width (m)

disp(['initial_area: ', num2str(initial_area)]);disp(['L_0: ', num2str(L_0)]);disp(['A_0: ', num2str(A_0)]);disp(['stretch_factor: ', num2str(stretch_factor)]);disp(['epsilon_axial: ', num2str(epsilon_axial)]);disp(['epsilon_transverse: ', num2str(epsilon_transverse)]);
% Update cross-sectional width
A_t = A_0 * (1 + epsilon_transverse); % Updated width (m)
disp(['A_t: ', num2str(A_t)]);
% Assuming unit thickness in z-direction
thickness = 1; % in meters

% Calculate expected force from stress-strain relationship
F_t = E * epsilon_axial * A_t * thickness; % Force in Newtons (N)
disp(['E: ', num2str(E), ' Pa']);disp(['epsilon_axial: ', num2str(epsilon_axial)]);disp(['A_t: ', num2str(A_t)]);disp(['thickness: ', num2str(thickness)]);disp(['F_t: ', num2str(F_t)]);
% Extract adhesion force from simulation
F_x = squeeze(F_store(time_index, :, 1)); % x-component of force at current time

% Total adhesion force in x-direction
F_adhesion = max(F_x); % Sum over all points (N)
disp(['max(F_x): ', num2str(max(F_x))]);
disp(['max(F_adhesion): ', num2str(F_adhesion)]);
F_adhesion = max(F_x); % Use the mean force as the adhesion force
disp(['F_adhesion: ', num2str(F_adhesion)]);
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
