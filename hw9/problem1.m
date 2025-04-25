%{
    Requires:
        - symbolic toolbox
%}

clear all;


%% Geometry
% OD of coil spring
OD = 31; % [mm]
% diameter of the music wire (not of overall coil spring)
d = 2.5; % [mm]
D = OD - d; % [mm]

% Total coils 
N_t = 14;
% See Table 10-1, plain and ground ends
% Number of active coils
N_a = N_t - 1;

%% Material Properties

% Music Wire
assert(d >= 0.1 & d <= 6.5, "Wire diameter falls outside range on table, manual intervention "...
    + "is required. See Table 10-4.");

% See Table 10-4
m = 0.145;
A = 2211; % [MPa * mm^m]

S_ut = A/(d^m); % [MPa]

% See Table 10-6
% Assume this is before set removed as there's no mention of the set being removed
allowable_torsional_stress = 0.45*S_ut; % [MPa]

assert(d/25.4 <= 0.125 & d/25.4 >= 0.064, "Wire diameter falls outside range on table, manual "...
    + "intervention is required.See Table 10-5."); % comparison is in [in]

G = 81.0*10^3; % [MPa]

%% Spring characteristics

% See Slide 5 of 23 in 330_S25_Lecture23_Springs.pdf
% [mm^4 * MPa * mm^-3] = [MPa*mm] = [N/mm^2 * mm] = [N/mm]
k = d^4*G / (8*D^3 * N_a); % [N/mm]

fprintf("k   = %.3f N/m\n", k*10^3);

assert(k*10^3 >= 1000 & k*10^3 <= 2000, "Spring constant falls outside range given for "...
    + "responses. Something has gone wrong."); % comparison is in [N/m]


% Use length to determine closing force
% Since n_s = 1, have yield strength upon

syms F; % [N]

% See Slide 4 of 23 in 330_S25_Lecture23_Springs.pdf
C = D/d;
K_B = (4*C + 2)/(4*C - 3);

% See Slide 4 of 23 in 330_S25_Lecture23_Springs.pdf
% [N*mm/mm^3] = [N/mm^2] = [MPa]
tau_max = K_B * 8*F*D / (pi*d^3); % [MPa]

% We're cutting costs with this one
% "Safety" factor
n_s = 1;

eq = tau_max*n_s == allowable_torsional_stress;
F = vpasolve(eq); % [N]

fprintf("F   = %.3f  N\n", F);

assert(F >= 100 & F <= 200, "Force value falls outside range given for "...
    + "responses. Something has gone wrong."); % comparison is in [N]

syms L_0;

% [N / (N/mm)]
delta_L = F/k; % [mm]

% See Table 10-1
L_s = d*N_t; % [mm]

L_0 = L_s + delta_L; % [mm]

fprintf("L_0 = %.3f    m\n", L_0*10^-3);

assert(L_0*10^-3 >= 0.1 & L_0*10^-3 <= 0.2, "Free length value falls outside range given for "...
    + "responses. Something has gone wrong."); % comparison is in [m]


% See Slide 7 of 23 in 330_S25_Lecture23_Springs.pdf
% Since we're told plain, ground end, and table specifies ends supported by flat surfaces must 
%   be squared and ground, assume it's either both ends pivoted or one end clamped and one end free 
% Assume both ends pivoted
alpha = 1;

if (L_0 < 2.75 * D / alpha)
    fprintf("N The spring is not predicted to buckle in service.\n");
else 
    fprintf("Y The spring is predicted to potentially buckle in service.\n");
end
