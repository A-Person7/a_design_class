% Requires symbolic toolbox

clear all;


% Material properties

Youngs_Modulus = 114 * 10^3; % [MPa]
Yield_Strength = 880; % [MPa]


% Overall Geometry
r_1 = 5; % [mm]
r_2 = 10; % [mm] 
L = 15; % [mm] 

% Loading conditions 
M_z = 20*10^3; % [N mm]
F_y = 900; % [N] 
F_x = 600; % [N]

% Find things as a function of x to be sure we have the worst case scenerio for von Mises 
%   equivalent stress
syms x; % [mm]
assume(x >= 0 & x <= L); % zero here has units of mm

% Utility function that takes in a symbolic expression a and substitutes in val for x, returning
%   a number to check calculations
ev = @(a, val) double(subs(a,x,val));

% Reactive forces
V = F_y; % [N]
N = - F_x; % [N]
M = M_z + F_y * (r_2 + x); % [N mm]


% Geometry of Cross Section at Cut
r = (r_1 * (L - x) + r_2 *x) / L; % [mm]

assert(abs(ev(r, 0) - r_1) < 0.001, "Check your radius expression.");
assert(abs(ev(r, L) - r_2) < 0.001, "Check your radius expression.");

A = pi * r^2; % [mm^2]
I = pi / 4 * r^4; % [mm^4]

% Bending case

% Minus sign on bending because we want to get the most stress, which is where bending
%   stress compounds with (compression/negative) axial stress
sigma_bending_case = -M * r / I + N / A; % [MPa]
tau_bending_case = 0; % [MPa]

% Shear case 
sigma_shear_case = N / A; % [MPa]
tau_shear_case = 4 * V / (3 * A); % [MPa]

% von Mises equivalent stress 
sigma_eq_bending_case = sqrt(sigma_bending_case^2 + 3 * tau_bending_case^2); % [MPa]
sigma_eq_shear_case = sqrt(sigma_shear_case^2 + 3 * tau_shear_case^2); % [MPa]


% Graphs

f = figure;
hold on;
x_num = [0:0.01:L]; % [MPa]
% Stress equivalent (von Mises stress)
seq_bending = double(subs(sigma_eq_bending_case, x, x_num)); % [MPa]
seq_shear = double(subs(sigma_eq_shear_case, x, x_num)); % [MPa]
plot(x_num, seq_bending, 'DisplayName', "Bending Case");
plot(x_num, seq_shear, 'DisplayName', "Shear Case");
legend
xlabel("x [mm]");
ylabel("von Mises stress [MPa]");
hold off;

% By inspection of graphs, bending is the worst case by far, and its equivalent stress strictly
%   decreases as x increases.

% Evaluate bending at worst case (x = 0 mm), cast to a number instead of a symbolic expression
sigma_eq_worst = double(subs(sigma_eq_bending_case, x, 0)); % [MPa]
fprintf("The worst stress experienced is %.3f MPa at the bottom of the narrowest region.\n",...
    sigma_eq_worst);
safety_factor = Yield_Strength / sigma_eq_worst;
fprintf("The factor of safety is %.3f.\n", safety_factor);

