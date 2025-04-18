% st is short for steel, ci is short for cast iron

clear all;



%% Material properties

E_st = 30; % [Mpsi]
E_ci = 14.5; % [Mpsi]
% Different steels still have roughly the same stiffness, assume bolt is steel
E_bolt = E_st; % [Mpsi]

%% Geometry

% Diameter of the bolt
d = 1/2; % [in]
% Diameter of head of the bolt, table lookup
D = 3/4; % [in]

% Length of the plates
l_st = 2; % [in]
l_ci = 1; % [in]

% Match notation used in class/slideshows
l = l_st + l_ci; % [in]

% See Table A-31, height of nut
H = 7/16; % [in]

% See Slide 26 of 51 in 330_S25_Lecture21_Fasteners_BasicsAndTensionJoints.pdf
L_min = l + H; % [in]

% assert(L_min > 3.4 && L_min < 3.5, "L_min value changed from original, manual intervention "...
%     + "required. Look up new L value in Table A-17 again and update this line.");
% 
% % Table A-17
% L = 3 + 1/2; % [in]

% Told to take ceiling of L to nearest 1/4 inch instead of using preferred sizes table.
L = 0.25*ceil(L_min*4); % [in]

fprintf("L   = %.3f in.\n", L);

assert(L <= 6, "L value changed, manual intervention required. Check Table 8-6.");

% Assume the minimum amount is threaded
% See Slide 51 of 51 in 330_S25_Lecture21_Fasteners_BasicsAndTensionJoints.pdf
L_t = 2*d + 1/4; % [in]
l_d = L - L_t; % [in]
l_t = l - l_d; % [in]

% Table 8-2
A_t = 0.1419; % [in^2]
A_d = pi*d^2/4; % [in^2]

%% Stiffness of Bolt

% [in^2 * in^2 * 10^6 * Mlbf / in^2 / (in^3)] = [Mlbf / in]
k_b = (A_d * A_t * E_bolt) / (A_d*l_t + A_t * l_d); % [Mlbf/in]

fprintf("k_b = %.3f Mlbf/in.\n", k_b);

%% Stiffness of Members

% Have three frusta, one in the steel until the halfway point, one in the steel that decreases 
%   until it hits the cast iron, and the other in the cast iron.
%   Let these be k_1, k_2, and k_3 respectively.

% Units of l/2 are [in], D are [in], units of d are [in], units of E_st are [Mlbf/in^2]
% Therefore, output units are guaranteed to be [Mlbf/in^2 * in] = [Mlbf/in]
k_1 = get_k_frustum(E_st, l/2, D, d); % [Mlbf / in]
% Using alpha = 30 (already assumed in function documentation)
k_2 = get_k_frustum(E_st, l/2 - l_ci, D + (l - l_ci)*tand(30), d); % [Mlbf / in]
k_3 = get_k_frustum(E_ci, l_ci, D, d); % [Mlbf / in]

% Stiffness sums inversely
k_m = (1/k_1 + 1/k_2 + 1/k_3)^(-1); % [Mlbf / in]

fprintf("k_m = %.3f Mlbf/in.\n", k_m);

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
