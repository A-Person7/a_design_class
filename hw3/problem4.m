clear all;

% All lengths in mm 

% See Figure A-15-9 in Shigley
d = 0.9*10; % [mm] 
D = 1.2*10; % [mm] 
c = d/2;
% Second area moment of inertia
I = pi*d^4 / 64;

% Shoulder fillet,
r = 3; % [mm]

% Moment applied
M = 10*10^3; % [N mm]

sigma_0 = M*c / I;

% Print out ratios to use to read from graph
fprintf("D/d is %.3f.\n", D/d)
fprintf("r/d is %.3f.\n", r/d)

K = 1.3;

sigma_max = K * sigma_0
