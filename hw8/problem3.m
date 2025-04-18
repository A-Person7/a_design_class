% Requires:
%   - symbolic toolbox

clear all;

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




%% Bolt Pattern

syms F; % [N]

% Worst case scenario for bolt is bolt A as shear term adds to moment term 

M = F*(r+c+R); % [N mm]
% Bolt A takes half of the moment about the center (M) 
shear_A = (F + M/2 * r)/A_d; % [MPa]
