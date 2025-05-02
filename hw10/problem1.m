%{
    Requires:
        - symbolic toolbox
%} 


clear all;

% script/cursive L on slideshow
script_L_D = 25; % [kh]
% rotational speed
n_D = 350; % [rev/min]
% See Slide 14 of 26 in 330_S25_Lecture24_Bearings.pdf
% [min/h * kh * rev/min] = [kilorev]
L_D = 60*script_L_D*n_D; % [kilorev]

radial_load = 2.5; % [kN]
application_factor = 1.2;

% Design reliability
R_D = 0.9;

% See Table 11-6, assume manufacturer 2 values (see Slide 11 of 26 in 330_S25_Lecture24_Bearings.pdf
%   for why)
L_10 = 1e3; % [kilorev]
x_0 = 0.02;
theta = 4.459;
b = 1.483;

% See Slide 14 of 26 in 330_S25_Lecture24_Bearings.pdf
% n_D is rot_speed
% multiple of life rating
% [kilorev/kilorev]
x_D = L_D / L_10;


print_and_validate("x_D", x_D, "", 500, 600);

% Ball bearing
a = 3;


% See Slide 9 of 26 in 330_S25_Lecture24_Bearings.pdf
% log is the natural log (not log_10). Run `help log` for more info
C_10 = application_factor * radial_load * ((x_D) / (x_0 + (theta-x_0)*(log(1/R_D))^(1/b)))^(1/a);

print_and_validate("C_10", C_10, "kN", 20, 30);

% See Table 11-2
assert(C_10 < 25.5 & C_10 > 19.5, "C_10 value changed, manual intervention is required for "...
    + "updating tabulated values. See Table 11-2.");

% Tabulated values, kept together to change them together
bore = 35; % [mm]
C_10 = 25.5; % [kN]

% fprintf("Bearing = 02-%2i mm\n", bore);
print_and_validate("Bearing", "02-"+bore, "mm", NaN, NaN);

% Update x_B
syms x_B;

% See Eq 11-2 in its modified form with x
% Substitute F_B with C_10 (equivalent), F_D by F*a_f
eq = C_10*(x_B)^(1/a) == radial_load*application_factor*(x_D)^(1/a);

% It's probably well posed
x_B = vpasolve(eq);

R = exp(-((x_B - x_0)/(theta - x_0))^(b));

print_and_validate("R", R*100, "%", 90, 95);

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
%       - Must be 7 characters or less
% - val -- the value of the variable in question
%       - can be a double, string, or character (array)
%       - will disable validation if it's of type string or char
% - units -- a string or character array representing the name of the units val is in 
%       - Must be 3 characters or less
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
        fprintf("%-7s = %-7s %-3s\n", name, val, units);
        return;
    end
 
    fprintf("%-7s = %-7.3f %-3s\n", name, val, units);
    
    % Check if range is NaN, if not, proceed with range checking
    if (~sum(isnan([range_low range_high])))
        assert(val >= range_low && val <= range_high, name + " is outside valid range given in"...
            + " problem statement. Manual intervention is required. Good luck.");
    end
end
