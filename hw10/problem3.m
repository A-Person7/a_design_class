clear all;

%% Parameters

diametral_pitch = 6; % [teeth/in]
% All full depth
N_driver = 22; % [teeth]
pressure_angle = 20; % [deg]
n = 1200; % [rev/min]
horsepower = 15; % [Hp], as one might expect
N_driven = 60; % [teeth]
face_width = 2; % [in]

%% Calculations

% Sanity check: [teeth / (teeth/in)] = [in]
pitch_dia = N_driver / diametral_pitch; % [in]

% [rev/min * rad/rev] = [rad/min]
omega = n * (2*pi); % [rad/min]

ft_per_in = 1/12; % [ft/in]

% See Slide 19 of 34 in 330_S25_Lecture25_Gears.pdf
V = pitch_dia/2 * omega * ft_per_in; % [ft/min]

% See Slide 19 of 34 in 330_S25_Lecture25_Gears.pdf
W_t = 33000 * (horsepower/V); % [lbf]

% For a spur gear, transverse diametral pitch is just teh diametral_pitch
% See Slide 24 of 34 in 330_S25_Lecture25_Gears.pdf
P_d = diametral_pitch; % [teeth/in]

% Engineering ahh approximation
K_o = 1;
K_s = 1;
K_m = 1;
K_B = 1;

% See Slide 26 of 34 in 330_S25_Lecture25_Gears.pdf
% For standard steel, cut or milled profile, full teeth, imperial units 
K_v = (1200 + V) / 1200;

% Use abs and an epsilon value to avoid floating point comparisons errors
assert(abs(pressure_angle - 20) < 0.00001 && abs(N_driver - 22) < 0.00001, "Tabulated values "...
    + "for Y will have to be re-read as gear profile has changed. Manual intervention required."...
    + " See Table 14-2.");

Y = 0.331;

% See Slide 25 of 34 in 330_S25_Lecture25_Gears.pdf
J = Y; 

psi_to_kpsi = 1e-3; % [kpsi/psi]

% See Slide 25 of 34 in 330_S25_Lecture25_Gears.pdf
% [lbf * teeth/in / in * kpsi/psi] = [lbf/in^2 * kpsi/psi] = [psi * kpsi/psi] = [kpsi]
sigma = W_t * K_o * K_v * K_s * P_d / face_width * K_m * K_B / J * 1e-3; % [kpsi]


%% Output Display
print_and_validate("σ", sigma, "kpsi", 6, 8);


%% Functions

% print_and_validate
% 
% Intended to print outputs in an aligned fashion directly applicable to the homework answer.
% Intended to be aesthetically pleasing, ended up kind of wonky, but oh well.
%
% If given a proper range (i.e. when Dr. Sargent includes bounds on a homework answer), this 
%   method can be used to validate and provide an obvious warning if values given are outside 
%   of this range. This behavior can be easily disabled.
% 
%
% Inputs:
% - name -- a string or character array containing the name of the variable.
%       - Must be 1 characters or less
% - val -- the value of the variable in question
%       - can be a double, any integer type, string, or character (array)
%       - will disable validation if it's of type string or char
% - units -- a string or character array representing the name of the units val is in 
%       - Must be 4 characters or less
% - range_low -- a double that contains the lower range value to check against 
%       - must match units of val
%       - set to NaN to disable all range checking (high and low)
%       - ignored if val is type string or char
%
% - range_high -- a double that contains the upper range value to check against 
%       - must match units of val
%       - set to NaN to disable all range checking (high and low)
%       - ignored if val is type string or char
%
% Side Effects:
%   - prints out formatted text indicating the value (and units) of the given variable
%   - 
%
% Throws/Assertations:
%   - assert fails when val is outside range_low and range_high (inclusive), but only if
%       - neither range_low nor range_high are NaN
%       - val is of type double
%
function print_and_validate(name, val, units, range_low, range_high)
    % Run `help strcmp` or `doc strcmp` for more info
    if (strcmp(class(val), 'string') || strcmp(class(val), 'char'))
        fprintf("%-1s = %-3s %-4s\n", name, val, units);
        return;
    end
    
    if (contains(class(val), "int"))
        % done instead of %i to avoid leading zero's
        fprintf("%-1s = %-3s %-4s\n", name, string(val), units);
    else
        % format as double
        fprintf("%-1s = %-3.3f %-4s\n", name, val, units);
    end
    
    % Check if range is NaN, if not, proceed with range checking
    if (~sum(isnan([range_low range_high])))
        assert(val >= range_low && val <= range_high, name + " is outside valid range given in"...
            + " problem statement. Manual intervention is required. Good luck.");
    end
end
