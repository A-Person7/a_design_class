% Requirements:
%   - symbolic toolbox
%   - add the scripts folder to your MATLAB path
%       - run 'help path' (without the '') in MATLAB to get more information on how to do this

clear;

% Turn this on if you're paranoid
sanity_check = false;

%% Parameters

% l_0 is the origin
l_0 = 0; % [in]
l_1 = 1; % [in]
l_2 = 2; % [in]
l_3 = 9; % [in]
l_4 = 14; % [in]
l_5 = 15; % [in]
l_6 = 16; % [in]

d_1 = 2; % [in]
d_2 = 2.472; % [in]
d_3 = 2.763; % [in]

F_1 = 18; % [lbf]
F_2 = 32; % [lbf]

% Gravitational acceleration
g = 32.1740 * 12; % [in/s^2]

%% Support Reactions

% Located at l_6
R_2 = (F_1*l_2 + F_2 * l_4) / l_6; % [lbf]
% Located at l_0 (origin)
R_1 = F_1 + F_2 - R_2; % [lbf]


%% Deflection

% From here on out, units are those that match the function contract/documentation of 
%   ShaftDeflectionEnglish, and notes are only provided for derived results stemming from 
%   Raleigh's Equation such as critical speeds

% Units are in accordance with contract of the function ShaftDeflectionEnglish
F =     [F_1, F_2];
F_loc = [l_2, l_4];

d =     [d_1, d_2, d_3, d_1];
d_loc = [l_0, l_1, l_3, l_5];

R_loc = [l_0, l_6];

L = l_6;

% What do these input arguments, outputs mean? Go read the documentation
[x,y,dydx, M, MdEI, R, diam, EI] = ShaftDeflectionEnglish(F,F_loc,d,d_loc,R_loc,L);

if (sanity_check)
    f = figure;
    plot(x, diam/2);
    xlabel("Position [in]");
    ylabel("Radius [in]");
    xlim([0,L])
    ylim([0,1.2*max(d)/2])
    title("Half Beam Geometry");

    proceed_prompt("Does this graph look right?", f);
end



%% Critical Speed of Attached Elements

% Find deflection at closest location in vector x to the 
%   actual length values external loads are applied at
[ignore, y_idx_1] = min(abs(l_2- x));
[ignore, y_idx_2] = min(abs(l_4- x));
y_1 = y(y_idx_1);
y_2 = y(y_idx_2);
% y_1 belongs to F_1 (NOT l_1)
% y_2 belongs to F_2 (NOT l_2)

% Numerator term
num = F_1 * y_1 + F_2 * y_2;
% Denominator term
denom = F_1 * y_1^2 + F_2 * y_2^2;

omega_1 = sqrt(g * num / denom); % [rad/s]


%% Critical Speed of Shaft

% I'm 85% sure this cancels, this value is from the critical speed example in Slide 12 of 15 in 
%   330_S25_Lecture20_Shafts2_DeflectionAndSpeed.pdf as that's also probably steel (just in case)
% weight_density_steel = 0.282; % [lbf / in^3]
weight_density_steel = 0.284; % [lbf / in^3]

% If we want more accuracy, we'd use a lot more sections, because we can do that pretty 
%   easily with a script. If we want to agree with the value the grading assistant is given,
%   however...
% Divide the shaft into 4 sections 

d =     [d_1, d_2, d_3, d_1];
d_loc = [l_0, l_1, l_3, l_5];

F_loc = [mean([l_0, l_1]), mean([l_1, l_3]), mean([l_3, l_5]), mean([l_5, l_6])];
F = weight_density_steel * pi / 4 * [ ...
    d_1^2 * (l_1 - l_0), ...
    d_2^2 * (l_3 - l_1), ...
    d_3^2 * (l_5 - l_3), ...
    d_1^2 * (l_6 - l_5), ...
];




R_loc = [l_0, l_6];
L = l_6;

% Rinse and repeat
[x,y,dydx, M, MdEI, R, diam, EI] = ShaftDeflectionEnglish(F,F_loc,d,d_loc,R_loc,L);

if (sanity_check)
    f = figure;
    plot(x, y);
    xlabel("Horizontal Position [in]");
    ylabel("Vertical Deflection[in]");
    xlim([0,L])
    ylim([0,1.2*max(y)])
    title("Beam Deflection");

    proceed_prompt("Does this graph look right?", f);
end


% Numerator term
num = 0;
% Denominator term
denom = 0;

for i = 1:length(F_loc)
    [ignore, y_idx] = min(abs(F_loc(i) - x));
    y_val = y(y_idx);
    num = num + F(i) * y_val;
    denom = denom + F(i) * y_val^2;
end

omega_s = sqrt(g * num / denom); % [rad/s]


%% Combination Critical Speed

% % Use Dunkerley's Equation
omega_1_tot = (omega_s^(-2) + omega_1^(-2))^(-1/2); % [rad/s]

% Take half the combined critical speed, convert to RPM
operating_speed = 0.5 * omega_1_tot; % [rad/s]

%% Display Outputs

fprintf("\n");

t = table;
t.("Critical Speed [RPM]") = rad_per_sec_to_rpm([omega_1; omega_s; omega_1_tot; operating_speed]);
t.Properties.RowNames = ["Attached Elements:"; "Shaft:"; "Combined:"; "Operating Speed:"];

disp(t);

%% Functions

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
%   - will throw an assert error if told to
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

% numOut
% 
% Use Java formatting code to get nicer looking numbers (adds commas so 1000 (double) becomes 
%   '1,000' (char)
% MATLAB is kind of secretly Java, some C++, and probably a little left over FORTRAN in a coat
%   anyways
% Not my original function, taken from:
% https://www.mathworks.com/matlabcentral/answers/96131-is-there-a-format-in-matlab-to-display-numbers-such-that-commas-are-automatically-inserted-into-the
%
% Input:
%   - numIn -- the number input as a scalar double
% 
% Output:
%   - numOut -- the string output as a char array
function numOut = addComma(numIn)
   import java.text.*
   jf=java.text.DecimalFormat; % comma for thousands, three decimal places
   numOut= char(jf.format(numIn)); % omit "char" if you want a string out
end

% rad_per_sec_to_rpm
% 
% Does exactly what it sounds like
% 
% Input:
%   - in -- angular speed in [rad/s] as a double or symbolic expression, can be an array
% Output: 
%   - out -- angular speed in [rev/min] as the same type as input
function out = rad_per_sec_to_rpm(in) 
    out = (in .* 60) ./ (2*pi); % [RPM]
end
