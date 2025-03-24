% Requires symbolic library
%   if you expect to have any figures open that you care about

clear all;

% Material Properties
compression_strength = 90 * 10^3; % [psi]
tensile_strength = 80 * 10^3; % [psi]


% Geometry
ID = 0.5; % [in]
OD = 0.75; % [in] 
L_x = 1; % [in]

% Use radius formula for circle, convert diameter to radius, and use linearity of integral 
%   to subtract center circle
I = pi/4 * (OD / 2)^4 - pi/4 * (ID/2)^4; % [in^4]
% Area
A = pi * ((OD/2)^2 - (ID/2)^2);


% Loading conditions
F = 250; % [lbf] 
M = 30; % [lbf in]


% Support reaction 
R_F = [0, F, 0]; % [lbf]
R_M = [0, M, F * L_x]; % [lbf in]

% angle of R_M with respect to yz plane, see diagram on paper for more information
phi = atan2d(R_M(2), R_M(3)); % [deg]

% Shear case 

% Take the integral, valid for points along the z axis (NOT IN GENERAL)
Q_shear_case = (OD^3 - ID^3) / 12; % [in^3]
% Valid for points along the z axis
b_shear_case = OD - ID; % [in]

V = dot(R_F, [0, 1, 0]);
tau_shear_case = V*Q_shear_case/(I * b_shear_case);
assert(tau_shear_case < 2*V / A & tau_shear_case > 4*V / (3*A), "Have a nonsensical tau value, check your Q calculation.");

sigma_z_shear_case = dot(R_M, [0, 1, 0]) * (OD/2) / I; % [psi]

% Function call is of the following form
% mohr_3d_fcn(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_xz)
[tau_max_shear_case, sigma_1_shear_case, sigma_2_shear_case, sigma_3_shear_case] = ...
    mohr_3d_fcn(0, 0, sigma_z_shear_case, 0, tau_shear_case, 0); % zeros are in psi

n_shear_case = (sigma_1_shear_case/tensile_strength - sigma_3_shear_case/compression_strength)^(-1);

% Bending case

sigma_z_bending_case = - (dot(R_M, [0, 1, 0]) * OD/2 * sind(phi)) / I; % [psi]
sigma_x_bending_case = - (dot(R_M, [0, 0, 1]) * OD/2 * cosd(phi)) / I; % [psi]

b_bending_case = OD*sind(phi);

syms y z
% We're returning to first principles with this one
% TODO -- check
Q_bending_case = double(int(2*int(y, z, 0, sqrt((OD/2)^2 - y^2)), y, OD/2*cosd(phi), OD/2)); % [mm]
tau_bending_case = V*Q_bending_case / (I * b_bending_case); % [MPa]
% What Dr. S (name partially redacted to reduce searchability) got, to check 
%   if this is a source of discrepancy
% tau_bending_case = 250; % [MPa]

% Function call is of the following form
% mohr_3d_fcn(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_xz)
% TODO -- fix
% [tau_max_bending_case, sigma_1_bending_case, sigma_2_bending_case, sigma_3_bending_case] = ...
%     mohr_3d_fcn(sigma_x_bending_case, 0, sigma_z_bending_case, 0, 0, tau_bending_case); % zeros are in psi

% Rotate axes to align with net resultant bendign moment to simplify calculations
[tau_max_bending_case, sigma_1_bending_case, sigma_2_bending_case, sigma_3_bending_case] = ...
    mohr_3d_fcn(0, sqrt(sigma_x_bending_case^2 + sigma_z_bending_case^2), 0, 0, tau_bending_case, 0); % zeros are in psi

n_bending_case = (sigma_1_bending_case/tensile_strength - sigma_3_bending_case/compression_strength)^(-1);

% I should really just bite the bullet and add a way of 'silently' calling mohr_3d_fcn without 
%   it printing everything and making a large graph. However, I don't want to retroactively 
%   edit submitted homework assignments, so, for now, this will have to do. 
% This just clears the console and closes the figure
clc;
close;

assert(n_bending_case < n_shear_case, "A horrible plight has fallen upon you, you've messed up "...
    + "this problem...");

fprintf("n = %.3f is our safety factor.\n", n_bending_case);
