clear all;

% NOTE: Shigley table numbers, figure numbers are with respect to the Tenth edition
%   This means there's a discrepancy with the slides where Table 6-4 on the slides (reliability
%   factor) is Table 6-5 in the Tenth Edition


% Notation:
%   _st is the AISI 1020 CD Steel
%   _al is the 319-T6 Aluminum

%% Geometry
b = 2.5*10; % [mm]
h = 0.5*10; % [mm]

%% Material Properties


% Table A-20 in Shigley
S_ut_st = 470; % [MPa]
% Table A-24
S_ut_al = 248; % [MPa]


%% Perfect Endurance limits/fatigue strength at 5E8 cycles
% Table A-24.b
S_f_5E8_prime_al = 69; % [MPa]
S_e_prime_st = NaN; % [MPa]

% Eq 6-10
if (S_ut_st <= 1400) % 1400 is in [MPa]
    S_e_prime_st = 0.5*S_ut_st;
elseif (S_ut_st > 1400) % 1400 is in [MPa]
    S_e_prime_st = 700; % [MPa]
else 
    assert(false, "You've messed something up, you're on your own on this one.");
end


fprintf("S_e_prime_st      = %6.3f MPa.\n", S_e_prime_st);
fprintf("S_f_5E8_prime_al  = %7.3f MPa.\n\n", S_f_5E8_prime_al);

% Now that we have our endurance limits/fatigue strength at 5E8, need Marin factors

%% k_a 

% Use Table 6-2

% Cold Drawn, in MPa
a_st = 3.04;
b_st = -0.217;

% Table A-24.b lists this aluminum alloy as cast (forged), in MPa
a_al = 54.9;
b_al = -0.758;

k_a_al = a_al*S_ut_al^b_al;
k_a_st = a_st*S_ut_st^b_st;

fprintf("k_a_st = %.3f \t k_a_al = %.3f\n", k_a_st, k_a_al);

%% k_b

% Same geometry/loading for both => same k_b;

% Axial loading
k_b_al = 1;
k_b_st = 1;

fprintf("k_b_st = %.3f \t k_b_al = %.3f\n", k_b_st, k_b_al);

%% k_c 

% Both are pure axial
k_c_al = 0.85;
k_c_st = 0.85;

fprintf("k_c_st = %.3f \t k_c_al = %.3f\n", k_c_st, k_c_al);

%% k_d 

% Operating temperature
T_op = 150; % [C]

% Use formula given in slides for temperature in Celcius
S_T_by_S_RT = 0.99 + 5.9*10^(-4)*T_op - 2.1 * 10^(-6)*T_op^2;

% Cap at max 1
S_T_by_S_RT = min(S_T_by_S_RT, 1);

% Explicitly told to use 0.95 for al
k_d_al = 0.95;
k_d_st = S_T_by_S_RT;

fprintf("k_d_st = %.3f \t k_d_al = %.3f\n", k_d_st, k_d_al);

%% k_e 

% We're told 97% reliability, but that's not listed in Table 6-4. Using my Engineering Intuition,
%   I'm just going to linearly interpolate
k_e_al = interp1([95, 99], [0.868, 0.814], 97);
k_e_st = k_e_al;

fprintf("k_e_st = %.3f \t k_e_al = %.3f\n", k_e_st, k_e_al);

%% k_f

% Metal spraying is present, see Slide 32 of 43 in 330_S25_Lecture14_Fatigue1_FullyReversed.pdf 
k_f_al = 0.86;
k_f_st = 0.86;

fprintf("k_f_st = %.3f \t k_f_al = %.3f\n", k_f_st, k_f_al);


%% Real Endurance limits/fatigue strength at 5E8 cycles

S_e_st = k_a_st * k_b_st * k_c_st * k_d_st * k_e_st * k_f_st * S_e_prime_st;
S_f_5E8_al = k_a_al * k_b_al * k_c_al * k_d_al * k_e_al * k_f_al * S_f_5E8_prime_al;

fprintf("\nS_e_st      = %6.3f MPa.\n", S_e_st);
fprintf("S_f_5E8_al  = %7.3f MPa.\n", S_f_5E8_al);
