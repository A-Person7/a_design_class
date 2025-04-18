%  ________________________________________
% / st is short for steel, ci is short for \
% \ cast iron                              /
%  ----------------------------------------

clear all;

%% Material Properties

E_st = 207 * 10^3; % [MPa]
E_ci = 100 * 10^3; % [MPa]

% Assume bolt is also steel, all steels have roughly the same stiffness
E_bolt = E_st; % [MPa]

% ISO 9.8, M12 x 1.75
% Proof strength, see Table 8-11
S_p = 650; % [MPa]

%% Geometry

% See diagram on problem statement
% Using _geom to differentiate from similarly named bolt terms
A_geom = 20;  % [mm]
B_geom = 20;  % [mm]
C_geom = 100; % [mm]
D_geom = 150; % [mm]
E_geom = 200; % [mm]
F_geom = 300; % [mm]

% Number of bolts
N = 10;

% Cap head area, defined by gasket
A_c = pi * (D_geom/2)^2; % [mm^2]

%% Loading

p_g = 6; % [MPa]

%% Bolt Geometry

d = 12; % [mm]
% See Table A-31
D = 18; % [mm]

% Height of nut, see Table A-31
H = 10.8; % [mm]
L_min = A_geom + B_geom + H; % [mm]

% See Table A-17
assert(L_min > 50 && L_min < 60, "L_min value changed from original, manual intervention "...
    + "required. Look up new L value in Table A-17 again and update this line.");
L = 60; % [mm]

fprintf("L   = %.3f mm\n", L);

% 1.75mm pitch is coarse, see Table 8-1
A_t = 84.3; % [mm^2]
A_d = pi*d^2/4; % [mm^2]

l = A_geom + B_geom; % [mm]

assert(L <= 125 & d <= 48, "Bolt changed, manual intervention required. Check slide 51 "...
    + "of 51 in 330_S25_Lecture21_Fasteners_BasicsAndTensionJoints.pdf.");

% Assume the minimum amount is threaded
% See Slide 51 of 51 in 330_S25_Lecture21_Fasteners_BasicsAndTensionJoints.pdf
L_t = 2*d + 6; % [mm]
l_d = L - L_t; % [mm]
l_t = l - l_d; % [mm]


%% Mechanics


k_b = (A_d * A_t * E_bolt) / (A_d*l_t + A_t * l_d); % [MN / mm]

assert(A_geom == B_geom, "Thickness of materials isn't the same, manual intervention is required "... 
    + "to update frustums.");

% Output units are guaranteed to be [MPa * mm] = [N/mm^2 * mm] = [N/mm]
k_st = get_k_frustum(E_st, A_geom, D, d); % [N/mm]
k_ci = get_k_frustum(E_ci, B_geom, D, d); % [N/mm]

% Stiffness sums inversely
k_m = (1/k_st + 1/k_ci)^(-1); % [N/mm]


C = k_b / (k_b + k_m);

% Proof load
F_p = A_t * S_p; % [N]
% Explicitly told non-permanent connection
F_i = 0.75*F_p; % [N]

% Assume no external pressure/given pressure is gauge pressure
P = p_g * A_c / N;

% See Slide 22 of 51 in 330_S25_Lecture21_Fasteners_BasicsAndTensionJoints.pdf
n_p = S_p * A_t / (C*P + F_i);

fprintf("n_p = %.3f.\n", n_p);

% See Slide 41 of 51 in 330_S25_Lecture21_Fasteners_BasicsAndTensionJoints.pdf
% [N/mm^2 * mm^2 - N / N] = [-]
n_L = (S_p * A_t - F_i) / (C * P);

fprintf("n_L = %.3f.\n", n_L);


% See Slide 41 of 51 in 330_S25_Lecture21_Fasteners_BasicsAndTensionJoints.pdf
% [N / N] = [-]
n_0 = F_i / (P * (1 - C) );

fprintf("n_0 = %.3f.\n", n_0);

%% Functions 


% get_k_frustum
% 
% Returns the stiffness of a frustum element within a homogenous material. Implicitly assumes
%   alpha = 30 degrees.
%
% See Figure 8-15, Equation 8-20 in Shigley's Mechanical Engineering Design, 10th Edition`
%
% Inputs:
%   - E -- Young's Modulus, scalar double, any units
%   - t -- the height of the frustum, scalar double
%       - must have the same units as d, D
%   - D -- the diameter of the top/bottom of the frustum (smaller of the two), scalar double
%       - must have the same units as t, d
%   - d -- the diameter of the hole cut from the center of the frustum, scalar double
%       - must have the same units as t, D
%
% Outputs:
%   k -- scalar, double, same units of E * units of d (e.g. lbf/in^2 * in = lbf/in)
%
function k = get_k_frustum(E, t, D, d)
    log_arg = ((1.155*t + D - d)*(D + d)) / ((1.155*t + D + d)*(D - d));

    % Recall that log in MATLAB is the natural log, not log_10 (see `help log` for more)
    % see Eq 8-20
    k = 0.5774*pi*E*d / log(log_arg);
end
