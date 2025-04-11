% Requires:
%   - symbolic toolbox

clear all;

% Set this to true if you're unsure of the results and want additional ways to validate things
% https://www.youtube.com/watch?v=BOTIIw76qiE
paranoid = false;


safety_factor = 2.5;


%% Geometry

% pitch diameter of the spur gear
pitch_dia = 150; % [mm]

% Length from B to C
center_length = 250; % [mm]
% Length from C to D
end_length = 100; % [mm]

force_angle = 20; % [deg]

%% Loading Conditions

T_A = 340 * 10^3; % [N mm]

F = T_A / (cosd(force_angle) * pitch_dia/2); % [N]


%% Material Properties
S_y = 420; % [MPa]
S_ut = 560; % [MPa]

%% Support Reactions

% Support reactions at the left and right bearings, respectively
syms R_B R_C; % [N]
support_rxn_eqs = [ ...
    R_B + R_C == F, ...
    R_C * center_length == F * (center_length + end_length) ...
];

sln = solve(support_rxn_eqs);
R_B = double(sln.R_B); % [N]
R_C = double(sln.R_C); % [N]

%% Shear and Bending Moments

% Let x be the distance, in mm, from B (0 mm) towards C
syms V(x) M(x) % [N], [N mm], respectively
% Dummy variable, standin for x
syms t; % [mm]
V(x) = R_B*heaviside(x) + R_C * heaviside(x - center_length) ...
    - F * heaviside(x - (center_length + end_length));

M(x) = int(subs(V(x), x, t), t, 0, x);

if (paranoid) 
    % For display
    extra_length = 10; % [mm]
    x_num = -extra_length:0.1:(center_length + end_length + extra_length); % [mm]
    v_num = double(subs(V, x, x_num)); % [N]
    m_num = double(subs(M, x, x_num)); % [N mm]

    f = figure;
    subplot(2, 1, 1);
    plot(x_num, v_num)
    xlabel("Position [mm]");
    ylabel("Shear [N]");
    title("Shear Diagram");

    subplot(2, 1, 2);
    plot(x_num, m_num)
    xlabel("Position [mm]");
    ylabel("Moment [N mm]");
    title("Bending Moment Diagram");

    proceed_prompt("Does this graph look right?", f);
end

%% Part a


% The diameter of the center length of shaft that we're interested in
syms dia; % [mm]
assume(dia, "positive");
assume(dia > 25 & dia < 75); % See answer range, both numbers have units of [mm]

% Ignore stress concentrations here because, uh, hmmm, erm...
%   Steel is ductile enough and this is static failure

% By inspection, worst case is clearly at x = 250mm
% It has 
%   - the largest momement (magnitude)
%   - the largest shear (magnitude) 
%   - stress concentration 

x_very_bad = 250; % [mm]

V_max = double(subs(V, x, x_very_bad)); % [N]
M_max = double(subs(M, x, x_very_bad)); % [N mm]

% Geometric properties of Cross Section
% Polar Moment of Inertia
polar_moment_of_inertia = pi * dia^4 / 32; % [mm^4]
% Use the linearity of integrals to our advantage
moment_of_inertia = 0.5 * polar_moment_of_inertia; % [mm^4]
cross_sec_area = pi * (dia/2)^2; % [mm^2]

% Torsional stress component is the same for both 
tau_torsion = T_A * (dia/2) / polar_moment_of_inertia; % [MPa]


% Now have two cases, max shear and max moment
% Max Shear
tau_shear = 4*V_max / (3*cross_sec_area); % [MPa]
tau_combined = tau_torsion + tau_shear; % [MPa]

% Dealing with plane stress, use equation in Slide 16 of 32 in 330_S25_Lecture13_StaticFailure.pdf
von_mises_shear_case = sqrt(3*tau_combined^2); % [MPa]

% Max Moment 
sigma = M_max * (dia/2) / moment_of_inertia; % [MPa]

% Dealing with plane stress, use equation in Slide 16 of 32 in 330_S25_Lecture13_StaticFailure.pdf
von_mises_moment_case = sqrt(sigma^2 + 3*tau_torsion^2); % [MPa]

% Solve
d_max_moment = double(solve(S_y / von_mises_moment_case == safety_factor)); % [mm]
d_max_shear = double(solve(S_y / von_mises_shear_case == safety_factor)); % [mm]

if (d_max_shear > d_max_moment)
    fprintf("For static failure, the maximum shear case governs.\n");
elseif (d_max_moment > d_max_shear) 
    fprintf("For static failure, the maximum moment case governs.\n");
else 
    % No way two different double precision numbers are going to be EXACTLY the same...
    fprintf("???");
end

d_static = max(d_max_shear, d_max_moment); % [mm]

fprintf("The minimum allowable diameter based on static yielding is %.2f mm.\n", d_static);

fprintf("\n");

%% Part b
% Analyze stress concentrations at the bearing shoulders or at the gear
% From Table 7-1, for first iteration
K_t = 2.7;
get_diam_from_fatigue(K_t, S_ut, tau_torsion, sigma, safety_factor);

% Use Figure A-15-9 to iterate
% This value was grabbed manually, needs to be updated manually if anything changes
K_t = 2.3;
get_diam_from_fatigue(K_t, S_ut, tau_torsion, sigma, safety_factor);



%% Functions

% get_diam_from_fatigue
% 
% Gets the minimal acceptable diameter of the shaft due to fatigue, assuming worst case scenario 
%   for stress is at max bending stress (not shear stress).
% 
% This function uses metric units, requires the symbolic toolbox, and requires the caller to expect 
%   the diameter to converge between a range of 25 mm and 75 mm. It is the caller's responsibility
%   to ensure the setup is such that the output converges, this function failing to converge is 
%   not well documented and can lead to unexpected behavior.
%
% Inputs:
%   K_t -- the current stress concentration factor to use for iteration as a scalar double 
%       - See Table 7-1 and Figure A-15-9
%   S_ut -- the ultimate strength of the material in question, in [MPa] as a scalar double 
%   tau_torsion -- the shear due to torsion, in [MPa], as a scalar symbolic expression function 
%       of dia (which is the diameter in [mm])
%   sigma -- the axial stress due to bending, in [MPa], as a scalar symbolic expression function 
%       of dia (which is the diameter in [mm])
%   safety_factor -- the design fator of safety as a scalar double
%
% Outputs:
%   - d_fatigue -- the minimum acceptable diameter based on the given K_t as a scalar double in [mm]
%
% Side Effects:
%   - prints output to the console
function d_fatigue = get_diam_from_fatigue(K_t, S_ut, tau_torsion, sigma, safety_factor)
    syms dia; % [mm]
    assume(dia, "positive");
    assume(dia > 25 & dia < 75); % See answer range, both numbers have units of [mm]
    
    % Assume that we want infinite life
    % Use moment case as worst case scenario, based on last part.
    
    % Ignore stress concentrations in the mean stress as it's a ductile material.
    sigma_m = sqrt(3*tau_torsion^2); % [MPa]
    
    % Not a particularly insightful formula...
    sigma_a_0 = sqrt(sigma^2); % [MPa]
    
    % To avoid magic numbers floating around
    mm_to_in = 1/25.4; % [in/mm]
    MPa_to_ksi =  0.1450377377; % [ksi/MPa]
    
    
    % Assume r/d = 0.02
    % Radius of fillet, in inches to match curve fit eq
    r_in = dia * 0.02 / mm_to_in; % [in]
    
    sqrt_a = 0.246-3.08e-3*(S_ut*MPa_to_ksi)+1.51e-5*(S_ut*MPa_to_ksi)^2- ...
        2.67e-8*(S_ut*MPa_to_ksi)^3; % [in^0.5]
    
    % I could use the graph, or I could use the Eq in Slide 28 of 43 in 
    %   330_S25_Lecture14_Fatigue1_FullyReversed.pdf
    q = 1 / (1 + sqrt_a / sqrt(r_in));
    
    K_f = 1 + q * (K_t - 1);
    
    sigma_a = K_f * sigma_a_0; % [MPa]
    
    endurance_strength = get_endurance_strength(S_ut, dia); % [MPa]
    
    % Since sigma_m >= 0 MPa, get this equation 
    % Don't even bother doing this symbolically, go directly to numerical solver
    % Read the MATLAB docs
    d_fatigue = vpasolve(sigma_a / endurance_strength + sigma_m / S_ut == 1/safety_factor); % [mm]
    
    fprintf("The minimum allowable diameter based on fatigue for a K_t value of " ...
        + "%.2f is %.2f mm.\n", K_t, d_fatigue);
end

% Gets the endurance strength of steel given the following conditions:
%   - as machined
%   - room temperature
%   - 50% reliability
%   - no miscellaneous effects
%
% The resultant value is only valid for comparison against von Mises stress
% In metric units.
%   
% ALL INPUTS, OUTPUTS MUST BE SCALARS
% 
% Inputs:
%   - S_ut - the ultimate (tensile) strength of the steel at room temperature [MPa]
%   - dia -- the diameter of the shaft, in [mm], as a double or as a symbolic expression
%       - must be between 7.62 mm and 254 mm exclusive for either (use asssumptions for symbolic)
%   
%  Outputs:
%   - endurance_strength -- the endurance strength (S_e) to compare against von Mises equivalent
%       stress [MPa]
%       - may be either symbolic expression (function of dia) or of type double depending on 
%           type of dia
function endurance_strength = get_endurance_strength(S_ut, dia)
    % See Eq 6-8
    S_e_prime = 0.5*S_ut; % [MPa]
    if (S_ut > 1400) % 1400 is in units of [MPa]
        S_e_prime = 700; % [MPa]
    end

    % See Table 6-2
    k_a = 3.04*S_ut^(-0.217);

    % See Slide 21 of 43 in 330_S25_Lecture14_Fatigue1_FullyReversed.pdf
    k_b = (dia/7.62)^(-0.107)*(heaviside(dia - 7.62) - heaviside(dia - 51))...
        + 1.51*dia^(-0.157)*(heaviside(dia - 51) - heaviside(dia - 254));
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

    endurance_strength = k_a * k_b * k_c * k_d * k_e * k_f * S_e_prime; % [MPa]
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

    % Never a bad time for a good old fashioned Linux message
    assert(~not_ok, "Bailing out, you are on your own. Good luck!");

    fprintf("Moving on...\n");
end
