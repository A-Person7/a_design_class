% Requires:
%   - symbolic toolbox to be installed


clear all;

%% Geometric Parameters

% Shaft diameters
d_1 = 1; % [in]
d_2 = 1.181; % [in]
d_3 = 1.7; % [in]
d_4 = 1.75; % [in]
d_5 = 2; % [in]
d_6 = 1.4; % [in]
d_7 = 1.181; % [in]

% Fillet radii
R_l = 1/16; % [in]
R_r = 1/32; % [in]
r_3 = 1/8; % [in]
r_4 = 0.1; % [in]
r_5 = 1/8; % [in]
r_6 = 1/32; % [in]


% Diameter change positions
l_0 = 0; % [in]
l_1 = 2; % [in]
l_2 = 2.75; % [in]
l_3 = 8.5; % [in]
l_7 = 12.87; % [in]
l_5 = l_7 - 2.2; % [in]
l_4 = l_5 - 0.485; % [in]
l_6 = l_7 - 0.75; % [in]

% Force locations 
% TODO -- CHECK. These values are not fully defined by the original problem statement, and are 
%   best guesses
L_1 = 1; % [in]
L_2 = 2.37; % [in]
L_3 = 9.35; % [in]
L_4 = l_7 - 0.375; % [in]


%% External Forces
F_radial = 230; % [lbf]
F_tangential = 633; % [lbf]


%% Solve Support Reactions

F = sqrt(F_radial^2 + F_tangential^2); % [lbf]

% Left (on diagram) bearing net reactive force, and right bearing reactive force, respectively
syms R_l R_r; % [lbf]

support_rxn_eqs = [ ...
    R_l + R_r == F, ... % sum forces in direction of net force)
    R_l*L_2 + R_r*L_4 == F*L_3 % sum moments from l_0 (origin)
];

sln = solve(support_rxn_eqs);
R_l = double(sln.R_l); % [lbf]
R_r = double(sln.R_r); % [lbf]
