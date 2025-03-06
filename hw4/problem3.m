% Requires Symbolic toolbox

clear;

% Inch, Pound, Second

% Since a lot of quantities don't change throughout the different cases and I'll need 
%   them for a few functions, put them all in a struct called param to make passing them 
%   as an argument easier. (Better than remembering to put globals everywhere)

% Safety factor
param.n_d = 2.5;
% If you're reading this, it means I either got the problem wrong, or you're just curious how I 
%   approached it. I'm hoping it's the later...
% Or you somehow stumbled on my GitHub account and were curious what other stuff I had

% Area of piston
A = pi * (2/2)^2; % [in^2]
% Compressive load
param.P = 1500 * A; % [lbf] = [psi * in^2]

% Table 4-2 in Shigley
param.C = 1;

% Young's Modulus, Table A5
param.E = 30 * 10^6; % [lbf/in^2] = [psi]
% Yield Strength, AISI 1030 HR Steel, Table A20
param.S_y = 37.5 * 10^3; % [psi]


[d_a, L_a] = get_min_distance_euler(50, param);
[d_b, L_b] = get_min_distance_euler(16, param);

fprintf("Assuming both are Euler buckling:\n");
fprintf("The minimum acceptable diameter for case a is %.3f inches.\n", d_a);
fprintf("The minimum acceptable diameter for case b is %.3f inches.\n", d_b);
fprintf("\n");

ratio_a = get_slenderness_ratio(d_a, L_a, param);
ratio_b = get_slenderness_ratio(d_b, L_b, param);

param.transition_slenderness_ratio = sqrt(2 * pi^2 * param.C * param.E / param.S_y);

if (ratio_a >= param.transition_slenderness_ratio)
    fprintf("Euler buckling criteria satisfied for case a.\n");
else 
    fprintf("Euler buckling criteria not satisfied for case a. Use Johnson buckling.\n");
end

if (ratio_b >= param.transition_slenderness_ratio)
    fprintf("Euler buckling criteria satisfied for case b.\n");
else 
    fprintf("Euler buckling criteria not satisfied for case b. Use Johnson buckling.\n");
end

fprintf("\nApplying Johnson buckling:\n");
[d_b, L_b] = get_min_distance_johnson(L_b, param);
fprintf("The minimum acceptable diameter for case b is %.3f inches.\n", d_b);

% Added manually, need to redo by hand if conditions change
fprintf("\nConverting to preferred sizes:\n");
% See Table A17 in Shigley
fprintf("The preferred size for case a is 1.20 inches.\n");
fprintf("The preferred size for case b is  3/4 inches.\n");



% Inputs:
%   - L -- the length of the rod, as drawn on paper [in]
%   - param -- the param struct generated earlier, all entries must be in base imperial units
% Outputs:
%   - D -- the diameter of the rod [in]
%   - L -- the length of the rod. Guaranteed to be the same as the input, done to
%       make it easier to keep together with the diameter
function [D, L] =  get_min_distance_euler(L, param)
    % Apply restriction d > 0 in to get physically possible output from solve
    syms d positive % [in]
    syms I % [in^4]
    I = pi/64 * (d)^4; % [in^4]
    P_cr = param.C * pi^2 * param.E * I / L^2; % [lbf]

    D = double(solve(P_cr/param.P == param.n_d)); % [in]
    L = L; % [in]
end

% Same contract as get_min_distance_euler
function [D, L] = get_min_distance_johnson(L, param)
    % Apply restriction d > 0 in to get physically possible output from solve
    syms d positive % [in]
    syms I % [in^4]
    syms A % [in^2]
    I = pi/64 * (d)^4; % [in^4]
    % Area
    A = pi * (d/2)^2; % [in^2]
    k = sqrt(I/A);
    P_cr = A * (param.S_y  - (param.S_y * L / (2*pi*k))^2 / (param.C * param.E));

    D = double(solve(P_cr/param.P == param.n_d)); % [in]
    L = L; % [in]
end

function l_by_k = get_slenderness_ratio(d, L, param)
    I = pi/4 * (d/4)^2; % [in^4]
    A = pi * (d/2)^2; % [in^2]
    k = sqrt(I/A); % [in]
    l_by_k = L / k;
end
