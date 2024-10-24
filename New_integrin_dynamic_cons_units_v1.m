function [integrin_free_new, integrin_bound_new] = New_integrin_cons_units_v2(...
    integrin_free, integrin_bound, F, k_0, k_1, k_m1, time_step_u, ...
    Concave_ind, F_0, G_f, D_f, D_b)

% Calculate the magnitude of the force at each point
absF = sqrt(sum(F.^2, 2));

% Update bound integrin concentration
delta_integrin_bound = ((k_0 + k_1 .* absF) .* integrin_free ...
    - k_m1 .* exp(-absF ./ F_0) .* integrin_bound ...
    - D_b .* integrin_bound) .* time_step_u;
integrin_bound_new = integrin_bound + delta_integrin_bound;

% Update free integrin concentration
delta_integrin_free = (G_f - D_f .* integrin_free ...
    - (k_0 + k_1 .* absF) .* integrin_free ...
    + k_m1 .* exp(-absF ./ F_0) .* integrin_bound) .* time_step_u;
integrin_free_new = integrin_free + delta_integrin_free;

% Ensure integrin concentrations are zero where there is no ECM
integrin_free_new = integrin_free_new .* Concave_ind;
integrin_bound_new = integrin_bound_new .* Concave_ind;

% Prevent negative concentrations
integrin_free_new(integrin_free_new < 0) = 0;
integrin_bound_new(integrin_bound_new < 0) = 0;
end
