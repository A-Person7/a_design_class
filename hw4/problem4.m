% Actually, this is Problem 2 again, but I figured it'd be more intuitive to call it Problem 4
%   Problem 4.132 in Shigley's 

% Requires Symbolic library

% Worst case scenario is when theta = 15 degrees
theta_worst = deg2rad(15); % [rad]

W = 300 * 9.80665; % [N]

P = W/4 ./ sin(theta_worst); % [N]

% Safety factor
n_d = 3.5;

% Let t be thickness
syms I(t) 

E = 207 * 10^3; % [MPa]
C = 1.4;
L = 350; % [mm]
w = 30; % [mm]
I = 1/12 * t^3 * w; % [mm^4, syms expression]

% Assume Euler buckling, check later
P_cr = C * pi^2 * E * I / L^2; %[N/mm^2 * mm^4 / mm^2 = N]

% Solve for t (only care about real root, cast to float)
t = double(solve(P_cr / P == n_d, 'Real', true)); % [mm]

fprintf("The minimum thickness is %.3f mm.\n", t);

% See Table A-17 in Shigley
t = max(t, 6.0); % Still [mm]

fprintf("The preferred thickness is %.1f mm.\n", t);

% See Table A-20
S_y = 180; % [MPa]

% Check if it's ok to use Euler buckling

transition_slenderness_ratio = sqrt(2 * pi^2 * C * E / S_y);

k = sqrt(double(subs(I, t))/(w*t)); % [mm]

% Compare
if (L/k > transition_slenderness_ratio)
    fprintf("Euler buckling assumptions satisfied.\n");
else 
    fprintf("Bad things have occurred.\n");
end

fprintf("Actual safety factor is %.3f.\n", double(subs(P_cr, t)) / P);

% Given that t < w and t is cubed in the expression, while C isn't, I'm just going to assume out 
%   of plane buckling governs (it does)
