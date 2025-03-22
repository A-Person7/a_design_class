% Requires symbolic library

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
tau_shear_case = V*Q/(I * b);
assert(tau_shear_case < 2*V / A & tau_shear_case > 4*V / (3*A), "Have a nonsensical tau value, check your Q calculation.");

sigma_z_shear_case = dot(R_M, [0, 1, 0]) * (OD/2) / I; % [psi]

% Function call is of the following form
% mohr_3d_fcn(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_xz)
[tau_max_shear_case, sigma_1_shear_case, sigma_2_shear_case, sigma_3_shear_case] = ...
    mohr_3d_fcn(0, 0, sigma_z_shear_case, 0, tau_shear_case, 0); % zeros are in psi

n_shear_case = (sigma_1_shear_case/tensile_strength - sigma_3_shear_case/compression_strength)^(-1);

% Bending case

sigma_z_bending_case = - (dot(R_M, [0, 1, 0]) * OD/2 * sin(phi)) / I; % [psi]
sigma_y_bending_case = - (dot(R_M, [0, 0, 1]) * OD/2 * cos(phi)) / I; % [psi]

b_bending_case = OD*sin(phi);
