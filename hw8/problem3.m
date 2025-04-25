% Requires:
%   - symbolic toolbox

clear all;

safety_factor = 2;

%% Geometry

% See paper diagram
t = 12;  % [mm]
h = 50;  % [mm]
r = 50;  % [mm]
c = 26;  % [mm]
R = 125; % [mm]

d = 10; % [mm]

A_d = pi*(d/2)^2; % [mm^2]

%% Material Properties

% ISO 5.8, M10 x 1.5
% Yield strength, see Table 8-11
S_y_bolt = 420; % [MPa]
% Use von Mises as it's what's used in Slide 10 of 18 in 
%   330_S25_Lecture22_Fasteners2_Shear&BoltPatterns.pdf
S_shear_bolt = 0.577*S_y_bolt; % [MPa]

%% Bolt Pattern

syms F; % [N]

% Worst case scenario for bolt is bolt A as shear term adds to moment term 

% Moment about the center
M = F*(r+c+R); % [N mm]

% Bolt A takes half of the moment about the center 
shear_A = (F + M/2 * r)/A_d; % [MPa]

eq_1 = shear_A == S_shear_bolt / safety_factor;

F_1 = double(solve(eq_1)); % [N]
