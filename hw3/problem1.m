% Parameters

clear all;

% Radius of the arm [m]
r = 0.75 * 10^-2;
% length of the arm [m]
L = 5 * 10^-2;
% 2nd Area Moment of Inertia of the arm
I = pi/4 * r^4; % [m^4]
J = 2*I; % [m^4]

F_1 = [0, 0, -500]; % [N]
r_F_1 = [2*L, 0, -L]; % [m]

F_2 = [400, 0, 0]; % [N]
r_F_2 = [2*L, 0, -L]; % [m]

M_1 = [-100, 0, 0]; % [Nm]
M_2 = [0, 0, 300]; % [N]

% Reactive force and moment at wall, respectively
R_F = -(F_1 + F_2); % [N]
R_M = -(cross(r_F_1, F_1) + cross(r_F_2, F_2) + ...
    M_1 + M_2); % [Nm]

% Only care about z direction of bending here
sigma_x = R_M(3) * r / I * 10^-6; % [MPa]
% No hoop stress
sigma_y = 0;  % [MPa]
% No "y" in -My/I here, therefore no z
sigma_z = 0;
tau_xy = 0;
tau_yz = 0;
% x direction of moment is torsion 
tau_xz = R_M(1)*r / J * 10^-6;;


mohr_3d_fcn(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_xz)
