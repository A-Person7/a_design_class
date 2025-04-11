% Requires:
%   - symbolic toolbox to be installed


% The smartest machine design sudent of today versus the hardest machine design problem of history
% https://www.youtube.com/watch?v=h_9RC8DCEPM

clear all;

graph_shear_bending = false;
graph_deflection = false;


%% Geometric Parameters
% See diagram on paper

% Shaft diameters
d_1 = 1; % [in]
d_2 = 1.181; % [in]
d_3 = 1.7; % [in]
d_4 = 1.75; % [in]
d_5 = 2; % [in]
d_6 = 1.4; % [in]
d_7 = 1.181; % [in]

% Fillet radii
r_1 = 1/16; % [in]
r_2 = 1/32; % [in]
r_3 = 1/8; % [in]
r_4 = 0.1; % [in]
r_5 = 1/8; % [in]
r_6 = 1/32; % [in]


% Diameter change positions
% l_0 to l_1 corresponds to d_1, l_1 to l_2 corresponds to d_2, etc.
l_0 = 0; % [in]
l_1 = 2; % [in]
l_2 = 2.75; % [in]
l_3 = 8.5; % [in]
l_7 = 12.87; % [in]
l_5 = l_7 - 2.2; % [in]
l_4 = l_5 - 0.485; % [in]
l_6 = l_7 - 0.75; % [in]

% Force locations (assume forces act in center of shoulder for their components)
L_1 = 1; % [in], fan keyway
L_2 = 2.375; % [in], left bearing
L_3 = mean([l_3, l_4 - 2*r_4]); % [in], gear keyway
L_4 = l_7 - 0.375; % [in], right bearing

% Driving gear pitch diameter
pitch_dia = 8; % [in]


%% External Forces
F_radial = 230; % [lbf]
F_tangential = 633; % [lbf]

% Torsion experienced within the shaft constantly between the keyways (L_1 and L_3)
torsion = F_tangential * pitch_dia/2; % [lbf in]


%% Solve Support Reactions

F = sqrt(F_radial^2 + F_tangential^2); % [lbf]

% Left (on diagram) bearing net reactive force, and right bearing reactive force, respectively
syms R_l R_r; % [lbf]

support_rxn_eqs = [ ...
    R_l + R_r == F, ... % sum forces in direction of net force)
    R_l*L_2 + R_r*L_4 == F*L_3 % sum moments from l_0 (origin)
];

sln = solve(support_rxn_eqs);

% These values are signed positive if directly opposing (anti-parallel) to F
R_l = double(sln.R_l); % [lbf]
R_r = double(sln.R_r); % [lbf]

fprintf("Left  bearing support reaction force is %.3f lbf.\n", R_l);
fprintf("Right bearing support reaction force is %.3f lbf.\n", R_r);

%% Internal Reactions

% Distance from the origin
syms x; % [in]
assume(x > L_1 & L_3 < 75); % Derived eq's only valid for 

syms V(x) M(x) % [lbf], [lbf in], respectively
% Dummy variable, standin for x
syms t; % [in]
V(x) = R_l*heaviside(x - L_2) + R_r * heaviside(x - L_4) - F * heaviside(x - L_3); % [lbf]
M(x) = int(subs(V(x), x, t), t, 0, x); % [lbf in]

if (graph_shear_bending)
    % Graph shear, bending moment diagram to validate expressions

    x_num = 0:0.001:l_7; % [in]
    v_num = double(subs(V, x, x_num)); % [lbf]
    m_num = double(subs(M, x, x_num)); % [lbf in]

    f = figure;
    subplot(2, 1, 1);
    plot(x_num, v_num)
    xlabel("Position [in]");
    ylabel("Shear [lbf]");
    title("Shear Diagram");

    subplot(2, 1, 2);
    plot(x_num, m_num)
    xlabel("Position [in]");
    ylabel("Moment [lbf in]");
    title("Bending Moment Diagram");

    fprintf("\n");
    proceed_prompt("Does this graph look right?", f);
end


%% Stresses, Fatigue


% AISI 1020 Cold Drawn Steel, see Table A-20
S_ut = 68; % [ksi]

% If it's not, you, yes, YOU, are responsible for modifying this and the relevant fatigue 
%   factors to change this script.
fprintf("\nAssuming bending is the worst case scenario.\n");


% Polar Moment of Inertia, dia is the diameter in [in]
polar_moment_of_inertia = @(dia) pi * dia.^4 ./ 32; % [in^4]
moment_of_inertia = @(dia) 1/2 * polar_moment_of_inertia(dia); % [in^4]

% Torsional stress, dia is in [mm]
tau_torsion = @(dia) torsion * (dia/2) / polar_moment_of_inertia(dia); % [psi]

% Ignore stress concentrations in the mean stress as it's a ductile material, dia is in [in]
% Valid for between L_1 and L_3, inclusive (points of interest for fatigue failure)
sigma_m = @(dia) sqrt(3*tau_torsion(dia)^2); % [psi]

% Technically, it's sqrt(sigma^2), but cut out the middle man
% x_val is a number for x, and is position measured from the origin [in]
sigma_a_0 = @(x_val, dia) abs(double(subs(M, x, x_val))  * (dia/2) ./ moment_of_inertia(dia)); % [psi]

% Apply stress concentration factors, assume worst case stress along keyways occurs in vicinity
%   of stress concentrations

fprintf("Assuming stress raiser location is worst case scenario.\n");

% See Figure A-15-9
% D/d = d_2 / d_1 = 1.181, r/d = r_1 / d_1 = 0.0625
sigma_a_at_1 = get_sigma_a(sigma_a_0(l_1, d_1), 1.75, r_1, S_ut); % [ksi]
endurance_strength_1 = get_endurance_strength(S_ut, d_1); % [ksi]
sigma_m_at_1 = sigma_m(d_1) * 10^-3; % [ksi]


incremental_x = 0.01; % [in]

% D/d = d_4 / (d_4 - 2*r_4) = 1.1290, r/d = r_4 / (d_4 - 2*r_4) = 0.0645
% Assume notch stress concentration factor is roughly that of a cylinder with the smaller diameter 
%   (conservative estimate)
% See Figure A-15-14
% Assume worst case scenario is in center of notch
sigma_a_at_3 = get_sigma_a(sigma_a_0(l_4-r_4, d_4 - 2*r_4), 2.1, r_4, S_ut); % [ksi]

endurance_strength_3 = get_endurance_strength(S_ut, d_4 - 2*r_4); % [ksi]
% Torsion is offloaded
sigma_m_at_3 = 0; % [ksi]

syms safety_factor;


% I could cut out the middle man and use vpasolve, but this ensures it' the proper answer 
%   and converges and all that
% I could also do the algebra, but...
safety_factor_at_1 = double(solve(...
    sigma_a_at_1 / endurance_strength_1 + sigma_m_at_1 / S_ut == 1/safety_factor)); 

safety_factor_at_3 = double(solve(...
    sigma_a_at_3 / endurance_strength_3 + sigma_m_at_3 / S_ut == 1/safety_factor)); 

fprintf("\n");
fprintf("The safety factor at the keyway in position %i is %.3f.\n", 1, safety_factor_at_1);
fprintf("The safety factor at the keyway in position %i is %.3f.\n", 3, safety_factor_at_3);
fprintf("\n");


%% Deflection Analysis

% Units are in accordance with contract of the function ShaftDeflectionEnglish
F =     [ F ];
F_loc = [L_3];

d =     [d_1 d_2 d_3 d_4 d_5 d_6 d_7];
d_loc = [l_0 l_1 l_2 l_3 l_4 l_5 l_6];

R_loc = [L_2, L_4];

L = l_7;

% Read the documentation for what these values mean, units of output
[x,y,dydx, M, MdEI, R, diam, EI] = ShaftDeflectionEnglish(F,F_loc,d,d_loc,R_loc,L);

if (graph_deflection)
    f = figure;
    subplot(3, 1, 1);

    plot(x, diam/2,'r');
    title("Half-Shaft Geometry");
    xlim([0,L]);
    ylim([0,1.2*max(d)/2]);
    ylabel("Radius");

    subplot(3,1,2);
    plot(x,y,'g');
    xlim([0,L]);
    title("Magnitude of Total Deflection");
    ylabel("Deflection [in]");
    xlabel("Position [in]");


    subplot(3,1,3);
    plot(x,dydx,'b');
    xlim([0,L]);
    title("Magnitude of total Slope")
    ylabel("Slope [rad]")
    xlabel("Position [in]");

    proceed_prompt("Does this graph look right?", f);
end

% Find deflection at closest location in vector x to the 
%   actual length values external loads are applied at
[ignore, x_idx_l_bearing] = min(abs(L_2- x));
[ignore, x_idx_r_bearing] = min(abs(L_4- x));
[ignore, x_idx_gear]      = min(abs(L_3- x));

dy_dx_l_bearing = dydx(x_idx_l_bearing);
dy_dx_r_bearing = dydx(x_idx_r_bearing);

dy_dx_gear = dydx(x_idx_gear);
y_gear = y(x_idx_gear);

% See Table 6-2
bearing_slope_limit = 0.003; % [rad]
gear_slope_limit = 0.0005; % [rad]
% Between 20 and 50 teeth/inch (worst case scenario)
gear_deflection_limit_fine = 0.003; % [in]
% Between 11 and 19 teeth/inch
gear_deflection_limit_middle = 0.005; % [in]
% Less than 10 teeth/inch
gear_deflection_limit_coarse = 0.010; % [in]

fprintf("\n");

fprintf("Slope value right bearing: %.3e rad.\n", dy_dx_r_bearing);
fprintf("\tAcceptable?                      "...
    + tern(abs(dy_dx_r_bearing) < bearing_slope_limit, "Yes", "No") + "\n");

fprintf("Slope value left bearing: %.3e rad.\n", dy_dx_l_bearing);
fprintf("\tAcceptable?                      "...
    + tern(abs(dy_dx_l_bearing) < bearing_slope_limit, "Yes", "No") + "\n");

fprintf("Slope value gear: %.3e rad.\n", dy_dx_gear);
fprintf("\tAcceptable? " + tern(abs(dy_dx_gear) < gear_slope_limit, "Yes", "No") + "\n");

fprintf("Deflection value gear: %.3e in.\n", y_gear);
fprintf("\tAcceptable worst case?           "...
    + tern(abs(y_gear) < gear_deflection_limit_fine, "Yes", "No") + "\n");
fprintf("\tAcceptable for 11-19 teeth/inch? "...
    + tern(abs(y_gear) < gear_deflection_limit_middle, "Yes", "No") + "\n");
fprintf("\tAcceptable for <10 teeth/inch?   "...
    + tern(abs(y_gear) < gear_deflection_limit_coarse, "Yes", "No") + "\n");


%% Functions 

% Gets the endurance strength of steel given the following conditions:
%   - Cold Drawn
%   - room temperature
%   - 50% reliability
%   - no miscellaneous effects
%
% The resultant value is only valid for comparison against von Mises stress
% In imperial units.
% 
% Uses the 11th Edition tabular values for k_a Marin factor (See Table 6-2 in 10th and 11th 
%   editions, there is a notable discrepancy).
%   
% ALL INPUTS, OUTPUTS MUST BE SCALARS
% 
% Inputs:
%   - S_ut - the ultimate (tensile) strength of the steel at room temperature [ksi]
%   - dia -- the diameter of the shaft, in [in], as a scalar double 
%   
%  Outputs:
%   - endurance_strength -- the endurance strength (S_e) to compare against von Mises equivalent
%       stress [ksi]
%       - may be either symbolic expression (function of dia) or of type double depending on 
%           type of dia
function endurance_strength = get_endurance_strength(S_ut, dia)
    % See Eq 6-8
    S_e_prime = 0.5*S_ut; % [ksi]
    if (S_ut > 200) % 200 is in units of [ksi]
        S_e_prime = 100; % [ksi]
    end

    % See Table 6-2
    % These values differ drastically between 10th and 11th Editions
    % 10th Edition:
    % k_a = 2.70*S_ut^(-0.265);
    % 11th Edition (also on slideshow)
    k_a = 2.00*S_ut^(-0.217);

    % See Slide 21 of 43 in 330_S25_Lecture14_Fatigue1_FullyReversed.pdf
    % It's faster to use a search and replace macro than properly convert to inches...
    k_b = ((dia*25.4)/7.62)^(-0.107)*(heaviside((dia*25.4) - 7.62) - heaviside((dia*25.4) - 51))...
        + 1.51*(dia*25.4)^(-0.157)*(heaviside((dia*25.4) - 51) - heaviside((dia*25.4) - 254));
    % Source on the above equation being right: trust me
    %   Also, I viewed it in Desmos to validate

    % Using von Mises stress, mentioned in the function documentation
    k_c = 1;
    % Room temperature, mentioned in function documentation
    k_d = 1;
    % 50% reliability, as mentioned in function documentation
    % See Table 6-4
    k_e = 1;
    % No miscellaneous factors
    k_f = 1;

    endurance_strength = k_a * k_b * k_c * k_d * k_e * k_f * S_e_prime; % [ksi]
end




% get_sigma_a
%
% Gets the alternating stress experienced in fatigue loading conditions due to stress 
%   concentration factors in units consistent with the endurance strength.
%
% Inputs:
%   - sigma_a_0 -- a double representing the pre stress concentration stress 
%       as a von Mises equivalent, [psi]
%   - K_t -- the K_t stress concentration factor to use as a double
%   - r -- the radius to use 
%   - S_ut -- the ultimate (tensile) strength of the material in question, in [ksi]
%
% Outputs:
%   - sigma_a -- the maximum alternating stress, in [ksi], experienced as a scalar double
function sigma_a = get_sigma_a(sigma_a_0, K_t, r, S_ut)
    sqrt_a = 0.246-3.08e-3*(S_ut)+1.51e-5*(S_ut)^2- ...
            2.67e-8*(S_ut)^3; % [in^0.5]
        
    % See Slide 28 of 43 in 330_S25_Lecture14_Fatigue1_FullyReversed.pdf
    q = 1 / (1 + sqrt_a / sqrt(r));
    
    K_f = 1 + q * (K_t - 1);

    sigma_a = sigma_a_0 * K_f * 10^-3;
end


% tern
% utility function to replicate a ternary operator
% Most proper languages have an operator where you can do condition ? expr1 : expr2 which is 
%   very nice for having quick, compact expressions
% https://en.wikipedia.org/wiki/Ternary_conditional_operator
% 
% Inputs:
%   - bool - a logical scalar 
%   - outA - first potential output, no restriction on type
%   - outB - second potential output, no restriction on type
% 
% Outputs:
% 
% Outputs:
%   - outA if bool is true, outB if bool is false
function out = tern(bool, outA, outB)
    if (bool)
        out = outA;
    else 
        out = outB;
    end
end



% Some code I wrote back in computer methods and modified slightly
%
% Inputs:
%   - query -- a string that represents a question to ask. y/n is automatically added after it
%   - varargin -- optional, an input variable that can be any number of function arguments 
%       - all must be figures or valid closable objects 
%       - these will all be closed by this function
% Outputs:
%   
% Side Effects:
%   - prints messages out into the console
%   - closes every varargin it gets 
% 
% Errors:
%   - will throw an assert error if told no
%
% prompts the user if they want to proceed
% inputs to this function are handles to be closed, and are an optional parameter
%   it's a variable argument, so you'd type it as proceed_prompt(f1, f2, f3);
% type 'y', 'Y', 'n', 'N', or hitting enter with nothing (treated as a continue)
function proceed_prompt(query, varargin) 
    % save user input as a proper string (instead of a char vector), and convert to lower case 
    %   to make matching easier
    str = lower(input(query + " y/n  ", "s"));

    % Is everything ok?
    % if the message is empty (you hit return), consider it as a yes, if n is present 
    %   treat it as a no
    not_ok = ~(contains(str, "y") || isempty(str)) || contains(str, "n");

    % varargin is a cell array instead of a vector, loop through all elements of it with this 
    %   loop
    for i = 1:length(varargin)
        f = varargin{i};
        % make sure f hasn't already been closed before attempting to close it
        if isvalid(f)
            close(f);
        end
    end

    assert(~not_ok, "Unfortunate.");

    fprintf("Moving on...\n");
end
