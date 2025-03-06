% See paper diagram for definition of x

clear;

% x is in meters, R is in Newtons, M is in Nm
syms M(x, R)

w = 5; % [N/m]
% syms w

% M = piecewise(0 <= x & x < 0.6, -w*x^2/2, ...
%     0.6 <= x & x < 1.5,... -w*x^2 / 2 + R*(x - 0.6), ...
%     1.5 <= x & x <= 2, -w*1.5*(x-0.75) + R*(x - 0.6));
% % To check this, you can copy the output of latex(M) (ignore the opening and closing '') and 
% %   paste it into a quick LaTeX -> png viewer (https://www.latex2png.com/)

% All Nm, x is in m
% M_a doesn't matter, but do anyways
M_a = -w*x^2 / 2; % 0 <= x < 0.6
M_b =  -w*x^2 / 2 + R*(x - 0.6); % 0.6 <= x < 1.5
M_c = -w*1.5*(x-0.75) + R*(x - 0.6); % 1.5 <= x < 2

% Can ignore 1/EI as it's set to 0 anyways. Bad practice, but can omit units as well 
int_a = int(M_a * diff(M_a, R), x);
int_b = int(M_b * diff(M_b, R), x);
int_c = int(M_c * diff(M_c, R), x);
indf_int = subs(int_a, x, 0.6) - subs(int_a, x, 0) + ...
    subs(int_b, x, 1.5) - subs(int_b, x, 0.6) + ...
    subs(int_c, x, 2.0) - subs(int_c, x, 1.5)

R = double(solve(indf_int == 0)); % [N]

fprintf("The reactive force at point C is %.3f Newtons.\n", R);


% integrand = M * diff(M, R); % [Nm^2]
% 
% % Indefinite integral
% integral = int(integrand, x); % [Nm^3]
% 
% % Evaluate at integral bounds (0m to 2m), but keep as symbolic
% integral_value = subs(integral, x, 2) - subs(integral, x, 0);  % [Nm^3]
% 
% % Solve for R, cast to a number
% % Since Nm^3 = N^2 m^3 / Units(R), we can see that R is in units of Newtons
% R = double(solve(integral_value == 0)); % [N]
% 
% fprintf("The reactive force at point C is %.3f Newtons.\n", R);
