clear all;

%% Geometry

% Total springs
N_b = 84;

d = 0.162; % [in]
OD = 1 + 1/2; % [in]
D = OD - d; % [in]

r_1 = D/2; % [in]
r_2 = 1/2; % [in]


%% Material Properties
% AISI 1065 OQ&T wire
% https://www.makeitfrom.com/compare/ASTM-A229-Oil-Tempered-Spring-Steel/SAE-AISI-1065-G10650-Carbon-Steel
% The above website says ASTM A229 and AISI 1065 Steel have nearly identical composition, therefore
%   use listed ASTM A229 values from Shigley to determine material properties.

assert(d >= 0.020 & d <= 0.50, "Wire diameter falls outside range on table, manual intervention "...
    + "is required. See Table 10-4."); % comparison is in [in]

% See Table 10-4
m = 0.187;
A = 147; % [ksi * in^m], where m is to the power of the variable m

S_ut = A/(d^m); % [ksi]

% Shigley does not include values for AISI 1065 or ASTM A229 Young's Modulus or Shear Modulus
% Instead, use listed values available online, which are very similiar to values for music wire 
%   and other alloys listed in Shigley's (as validation)
% https://www.azom.com/article.aspx?ArticleID=6575
E = 29007; % [ksi]
% When a pasta box says it cooks in 10-12 minutes, a wise person sets their timer for 11 minutes.
poissons_ratio = mean([0.27 0.3]);

G = E / (2*(1 + poissons_ratio)); % [ksi]

% Comparison is in Mpsi, check against similiar alloys in Shigley
assert(E*10^-3 > 27 & E*10^-3 < 30 & G*10^-3 > 11 & G*10^-3 < 12.5, "Your Young's Modulus and/or"...
    + " Shear Modulus values seem to deviate from what you'd expect of similiar alloys.");



%% Spring Characteristics

F_preload = 16; % [lbf]

C = D/d;
% See Slide 4 of 23 in 330_S25_Lecture23_Springs.pdf
K_B = (4*C+2)/(4*C-3);

% Eq 10-39
L_0 = 2*(D-d) + (N_b + 1)*d; % [in]

fprintf("L_0 = %6.3f in\n", L_0);

assert(L_0 >= 10 & L_0 <= 20, "Free length value falls outside range given for "...
    + "responses. Something has gone wrong."); % comparison is in [in]

% See Slide 12 of 23 in 330_S25_Lecture23_Springs.pdf
N_a = N_b + G/E;

% See Slide 4 of 23 in 330_S25_Lecture23_Springs.pdf
% [lbf * in * in^-3] = [psi]
tau_preload = K_B * 8*F_preload*D / (pi*d^3); % [psi]

fprintf("tau = %6.3f kpsi\n", tau_preload*10^-3);

assert(tau_preload*10^-3 >= 10 & tau_preload*10^-3 <= 20, "Free length value falls outside "...
    + "range given for responses. Something has gone wrong."); % comparison is in [ksi]

% So there's no magic numbers floating around (a bit late for that, I'll admit)
ksi_to_psi = 10^3;

% See Slide 5 of 23 in 330_S25_Lecture23_Springs.pdf
% [in^4 * lbf/in^2 * in^-3] = [lbf/in]
k = d^4 * G*ksi_to_psi / (8*D^3 * N_a); % [lbf/in]

fprintf("k   = %6.3f lbf/in\n", k);

assert(k >= 1 & k <= 5, "Spring constant estimate value falls outside range given for responses."...
    + " Something has gone wrong."); % comparison is in [lbf/in]

% Check body load for permanent deformation
% Assume no 'safety factor' (want most accurate estimate for permanent deformation that's the center 
%   of the distribution of failure)
% Assume static failure, use ends as they're the worst case
% See Slide 14 of 23 in 330_S25_Lecture23_Springs.pdf
S_allowable_torsion = 0.4*S_ut; % [ksi]
S_allowable_normal = 0.75*S_ut; % [ksi]

syms F; % [lbf]
% See Slide 4 of 23 in 330_S25_Lecture23_Springs.pdf
% [lbf * in * in^-3] = [psi]
tau = K_B * 8*F*D / (pi*d^3); % [psi]
F_max_body = vpasolve(tau == S_allowable_torsion * ksi_to_psi);


% Check hook load for permanent deformation

% See Slide 13 of 23 in 330_S25_Lecture23_Springs.pdf
C_1 = 2*r_1 / d;
% (K)_A, not 'just' K_A
K_A = (4*C_1^2 - C_1 -1) / (4*C_1*(C_1 -1));
% [lbf * in^-2]
sigma_A = F*(K_A * (16*D)/(pi*d^3) + 4/(pi*d^2)); % [psi]
C_2 = 2*r_2 / d;
% (K)_B, not 'just' K_B
K_B = (4*C_2 - 1) / (4*C_2 - 4);
% [lbf * in * in^-2] = [lbf/in^2]
tau_B = K_B * (8*F*D)/(pi*d^3); % [psi]

F_max_hook_normal = vpasolve(sigma_A == S_allowable_normal*ksi_to_psi); % [lbf]
F_max_hook_torsion = vpasolve(tau_B == S_allowable_torsion*ksi_to_psi); % [lbf]

% As of writing this, hook normal stress governs, but this picks the lowest maximum force before 
%   the sprign permanently deforms
F_max = min([F_max_body F_max_hook_torsion F_max_hook_normal]); % [lbf]

fprintf("F   = %6.3f lbf\n", F_max);

assert(F_max >= 80 & F_max <= 100, "Force to permanently deform spring value falls outside "...
    + "range given for responses. Something has gone wrong."); % comparison is in [lbf]

% [lbf * in/lbf]
deflection = (F_max - F_preload) / k; % [in]

fprintf("y   = %6.3f in\n", deflection);

assert(deflection >= 10 & deflection <= 20, "Deflection value falls outside "...
    + "range given for responses. Something has gone wrong."); % comparison is in [in]
