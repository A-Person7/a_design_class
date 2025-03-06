% Requires Symbolic toolbox

clear;

M = 0.3; % [Nm]
F = 20; % [N]

% Q_n is a dummy load (0N)
syms R_m_left R_m_right x Q_a Q_b % [Nm], [Nm], [m], [N], [N], respectively

% left has values for the left hand side (x \in [1m, 3m]) and right has values for x \in [0m, 1m)
%   where x is as-drawn on paper


E_left = 5*10^9; % [Pa]
h_left = 0.15; % [m]
d_left = 0.20; % [m]
I_left = 1/12 * h_left^3 * d_left; % [m^4]

E_right = 3*10^9; % [Pa]
h_right = 0.05; % [m]
d_right = 0.10; % [m]
I_right = 1/12 * h_right^3 * d_right; % [m^4]

% Valid for the left side of the beam, that 1 has units of meters
R_m_left = M + F*x + Q_b * x + Q_a * (x - 1); % [Nm]

% Valid for the right side of the beam
R_m_right = M + F*x + Q_b * x; % [Nm]

% Deflection at a:

% [Nm / (N/m^2 * m^4) * Nm/N] = [unitless] (checks out, multiplied by dt [m] to get deflection)
left_integrand_a = R_m_left / (E_left * I_left) * diff(R_m_left, Q_a);
left_integrand_a = subs(subs(left_integrand_a, Q_a, 0), Q_b, 0); % Swap out dummy loads with 0N
right_integrand_a = R_m_right / (E_right * I_right ) * diff(R_m_right, Q_a);
right_integrand_a = subs(subs(right_integrand_a, Q_a, 0), Q_b, 0); % Swap out dummy loads with 0N

left_integral_a = int(left_integrand_a); % [m]
right_integral_a = int(right_integrand_a); %[m]

% Evaluate both integrals at their bounds, cast to a float
delta_a = double(...
    subs(left_integral_a, x, 3) - subs(left_integral_a, x, 1) + ... % 3, 1 are in units of [m]
    subs(right_integral_a, x, 1) - subs(right_integral_a, x, 0) ... % 1, 0 are in units of [m]
); % [m]


% Conveniently, I can just copy and paste, and substitute _a with _b for deflection at b 
%   (except for substituting Q_a for 0N)
% Deflection at b:

% [Nm / (N/m^2 * m^4) * Nm/N] = [unitless] (checks out, multiplied by dt [m] to get deflection)
left_integrand_b = R_m_left / (E_left * I_left) * diff(R_m_left, Q_b);
left_integrand_b = subs(subs(left_integrand_b, Q_a, 0), Q_b, 0); % Swap out dummy loads with 0N
right_integrand_b = R_m_right / (E_right * I_right ) * diff(R_m_right, Q_b);
right_integrand_b = subs(subs(right_integrand_b, Q_a, 0), Q_b, 0); % Swap out dummy loads with 0N

left_integral_b = int(left_integrand_b); % [m]
right_integral_b = int(right_integrand_b); %[m]

% Evaluate both integrals at their bounds, cast to a float
delta_b = double(...
    subs(left_integral_b, x, 3) - subs(left_integral_b, x, 1) + ... % 3, 1 are in units of [m]
    subs(right_integral_b, x, 1) - subs(right_integral_b, x, 0) ... % 1, 0 are in units of [m]
); % [m]


% Since Q_a, Q_b are both signed such that +y direction is positive, so are their respective 
%   deflections.
fprintf("The deflection at a is %.3f mm in the positive y direction.\n", delta_a * 10^3);
fprintf("The deflection at b is %.3f mm in the positive y direction.\n", delta_b * 10^3);
