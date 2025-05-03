% Requires:
%   - symbolic toolbox

clear all;

safety_factor = 2;

%% Geometry

% See paper diagram
t = 12;  % [mm]
h = 50;  % [mm]
r = 50;  % [mm]
c = 26;  % [mm]
R = 125; % [mm]

% See Table A-7
bolt_bearing_length = 6.4; % [mm]

d = 10; % [mm]

A_d = pi*(d/2)^2; % [mm^2]

%% Material Properties

% ISO 5.8, M10 x 1.5
% Yield strength, see Table 8-11
S_y_bolt = 420; % [MPa]
% Use von Mises as it's what's used in Slide 10 of 18 in 
%   330_S25_Lecture22_Fasteners2_Shear&BoltPatterns.pdf
S_shear_bolt = 0.577*S_y_bolt; % [MPa]

% AISI 1006 steel
S_y_channel = 170; % [MPa]
% AISI 1015 steel
S_y_cantilever= 190; % [MPa]

%% Bolt Pattern

syms F; % [N]

% Worst case scenario for bolt is bolt A as shear term adds to moment term 

% Moment about the center
M = F*(r+c+R); % [N mm]

% Bolt A takes half of the moment about the center 
% tau_A = (F + M/2 * r)/A_d; % [MPa]
% Split the 'direct' shear amongst the three bolts equaly, take the moment
tau_A = (F/3 + M/(2*r))/A_d; % [MPa]

eq_1 = tau_A == S_shear_bolt / safety_factor;

F_1 = vpasolve(eq_1); % [N]


%% Bearing on Bolts

% Same force on the bolt (holding F equal) as last part of problem
eq_2 = S_y_bolt / safety_factor == (F/3 + M/(2*r)) / (d*bolt_bearing_length);

F_2 = vpasolve(eq_2); % [N]

%% Bearing on Channel

% Same as last part, except different material properties
eq_3 = S_y_channel / safety_factor == (F/3 + M/(2*r)) / (d*bolt_bearing_length);
F_3 = vpasolve(eq_3); % [N]


%% Bearing on Cantilever 

% deja vu
eq_4 = S_y_cantilever / safety_factor == (F/3 + M/(2*r)) / (d*t);
F_4 = vpasolve(eq_4); % [N]


%% Bending of Cantilever

% sigma = Mc/I

I = 1/12 * t * (r^3 - d^3); % [mm^4]

eq_5 = S_y_cantilever / safety_factor == F*(R+c) * (r/2) / I;
F_5 = vpasolve(eq_5); % [N]

%% Display 


F_governing = min([F_1 F_2 F_3 F_4 F_5]); % [N]

print_and_validate("Safe F", F_governing*1e-3, "kN", 1, 7.5);

print_and_validate("F_1", F_1*1e-3, "kN", 1, 7.5);
print_and_validate("F_2", F_2*1e-3, "kN", 1, 7.5);
print_and_validate("F_3", F_3*1e-3, "kN", 1, 7.5);
print_and_validate("F_4", F_4*1e-3, "kN", 1, 7.5);
print_and_validate("F_5", F_5*1e-3, "kN", 1, 7.5);


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
%       - Must be 6 characters or less
% - val -- the value of the variable in question
%       - can be a double, any integer type, string, or character (array)
%       - will disable validation if it's of type string or char
% - units -- a string or character array representing the name of the units val is in 
%       - Must be 5 characters or less
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
        fprintf("%-6s = %7s %-5s\n", name, val, units);
        return;
    end
    
    if (contains(class(val), "int"))
        % done instead of %i to avoid leading zero's
        fprintf("%-6s = %7s %-5s\n", name, string(val), units);
    else
        % format as double
        fprintf("%-6s = %7.3f %-5s\n", name, val, units);
    end
    
    % Check if range is NaN, if not, proceed with range checking
    if (~sum(isnan([range_low range_high])))
        assert(val >= range_low && val <= range_high, name + " is outside valid range given in"...
            + " problem statement. Manual intervention is required. Good luck.");
    end
end
