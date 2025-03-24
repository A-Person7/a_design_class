% Requires symbolic toolbox

clear all;


% Material Properties 
yield_strength = 542; % [MPa]


% Geometry
OD = 12.5; % [cm]
ID = 10.0; % [cm]
r_o = OD/2; % [cm]
r_i = ID/2; % [cm]

% Loading Conditions
internal_pressure = 10; % [MPa]


% Stress Calculations

% At worst case, r = r_i 
sigma_r = - internal_pressure * ( ((r_o/r_i)^2 - 1) / ((r_o / r_i)^2 - 1) ); % [MPa]
sigma_theta = internal_pressure * ( ((r_o/r_i)^2 + 1) / ((r_o / r_i)^2 - 1) ); % [MPa]
sigma_z = internal_pressure * ( (r_i)^2  / ( r_o^2 - r_i^2) ); % [MPa]

% No shear => these are our principle stresses

sigmas = sort([sigma_r, sigma_theta, sigma_z], 'descend'); % [MPa]
sigma_1 = sigmas(1); % [MPa]
sigma_2 = sigmas(2); % [MPa]
sigma_3 = sigmas(3); % [MPa]

tau_max = (sigma_1 - sigma_3) / 2; % [MPa]


% Max Shear 
n_max_shear = (yield_strength/2) / tau_max;


% Distortion Energy
sigma_eq = sqrt(0.5*( ...
    (sigma_1 - sigma_2)^2 + ...
    (sigma_2 - sigma_3)^2 + ...
    (sigma_3 - sigma_1)^2   ...
)); % [MPa]

n_distortion_energy = yield_strength / sigma_eq;


% Output
fprintf("Using Max Shear failure theory:         \tn = %.3f is our safety factor.\n", n_max_shear);
fprintf("Using Distortion Energy failure theory: \tn = %.3f is our safety factor.\n", n_distortion_energy);
