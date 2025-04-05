% Requires:
%   - symbolic toolbox

clear all;

% Set this to true if you're unsure of the results and want additional ways to validate things
% https://www.youtube.com/watch?v=BOTIIw76qiE
paranoid = true;


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
% V(x) = R_B*heaviside(x) - R_C * heaviside(x - center_length) + F * heaviside(x - (center_length - end_length));
V(x) = R_B*heaviside(x) + R_C * heaviside(x - center_length) ...
    - F * heaviside(x - (center_length + end_length));

M(x) = int(subs(V(x), x, t), t, 0, x);

if (paranoid) 
    % For display
    extra_length = 10; % [mm]
    x_num = -extra_length:0.1:(center_length + end_length + extra_length);
    v_num = double(subs(V, x, x_num));
    m_num = double(subs(M, x, x_num));

    f = figure;
    subplot(2, 1, 1);
    plot(x_num, v_num)
    xlabel("Position [in]");
    ylabel("Shear [N]");
    title("Shear Diagram");

    subplot(2, 1, 2);
    plot(x_num, m_num)
    xlabel("Position [in]");
    ylabel("Moment [N mm]");
    title("Bending Moment Diagram");

    proceed_prompt("Does this graph look right?", f);
end

%% Part a


% The diameter of the center length of shaft that we're interested in
syms dia; % [mm]


% Analyze stress concentrations at the bearing shoulders or at the gear

% From Table 7-1, for first iteration
K_t = 2.7;

% By inspection, worst case is clearly at x = 250mm
% It has 
%   - the largest momement (magnitude)
%   - the largest shear (magnitude) 
%   - stress concentration 

x_very_bad = 250; % [mm]

V_max = double(subs(V, x, x_very_bad)); % [N]
M_max = double(subs(M, x, x_very_bad)); % [N mm]



%% Part b

% Assume that we want infinite life




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
