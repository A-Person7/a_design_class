clear all;

%% Setup 

pressure_angle = 20; % [deg]
diametral_pitch = 6; % [teeth/in]

% Have 3 degrees of freedom,
%   Three of N_2, N_3, N_4, N_5 can be set free, with the remaining being whatever makes the 
%       proper gear ratio

% Through the magic of looking at the answer ranges, I can know when to stop my search
% Inclusive
% In SCREAMING_SNAKE_CASE as they're constants 
MIN_TEETH_NUM = 10; % [teeth]
MAX_TEETH_NUM = 225; % [teeth]
% Since it's a compound gear train, see Slide 14 of 34 in 330_S25_Lecture25_Gears.pdf
MAX_GEAR_RATIO = 100;

clearance = 0.5; % [in]
wall_thickness = 0.75; % [in]

% Includes N_2, N_3, N_4, N_5 in that order
gear_train = [NaN, NaN, NaN, NaN, NaN]; % [teeth]
gearbox_size = NaN; % [in]

GEAR_RATIO = 40;

% For floating point precision comparisons
% This very small number is used as a standin for zero
% See `help global` for more
global epsilon;
epsilon = 0.00001;

%% Brute Force

% Chose N_5 to be fixed, others to be free
% Brute force every potentially valid gear tooth combination 
% Might consider 10077696 combinations, of which only 7 are valid (maybe less)
for N_2 = MIN_TEETH_NUM:MAX_TEETH_NUM
    for N_3 = MIN_TEETH_NUM:MAX_TEETH_NUM
        for N_4 = MIN_TEETH_NUM:MAX_TEETH_NUM
            % Pick N_5 as that teeth value that forms the proper gear ratio 
            % As with many things gear or ratio related, this equation requires a 'magic 
            %   inverse' of the GEAR_RATIO
            N_5 = (1.0/GEAR_RATIO) * (N_2 * N_4) / N_3; % [teeth]

            % Since this is a loop, continue can be considered to 'prune' a gear train combination 
            %   if it isn't valid. See `help continue` for more

            % Prune too many teeth, too little teeth, and non-integer numbers of teeth, respectively
            if (N_5 > MAX_TEETH_NUM || N_5 < MIN_TEETH_NUM || ~is_int(N_5))
                continue;
            end

            % Quickly check if gear ratios are valid for compound spur gears
            if (N_4/N_5 > MAX_GEAR_RATIO || N_5/N_4 > MAX_GEAR_RATIO || ...
                    N_2/N_3 > MAX_GEAR_RATIO || N_3/N_2 > MAX_GEAR_RATIO)
                continue;
            end

            % Fail gears that'd need undercutting
            % full depth teeth
            k = 1;
            m_2_3 = max([N_2 N_3]) / min([N_2 N_3]);
            m_4_5 = max([N_4 N_5]) / min([N_2 N_3]);
    
            % I'm going to play fast and loose with driver vs. driven for a pair because this 
            %   loop setup (at this stage) could allow either of the two to be larger/smaller 
            % See Slide 10 of 34 in 330_S25_Lecture25_Gears.pdf
            N_min_2_3 = 2*k / ((1+2*m_2_3)*(sind(pressure_angle))^2) * ...
                (m_2_3 + sqrt(m_2_3^2 + (1+2*m_2_3)*(sind(pressure_angle)^2))); % [teeth]
            N_min_4_5 = 2*k / ((1+2*m_4_5)*(sind(pressure_angle))^2) * ...
                (m_4_5 + sqrt(m_4_5^2 + (1+2*m_4_5)*(sind(pressure_angle)^2))); % [teeth]

            if (N_2 < N_min_2_3 || N_3 <= N_min_2_3 || N_4 <= N_min_4_5 || N_5 <= N_min_4_5)
                continue;
            end

            % "War is Peace. Freedom is Slavery. Have the gear ratios of the sets be integer 
            %   values."
            %{ 
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@%.%%%%%%%%%%%%%%%%%%%%.@@@@@@@@@@.%%%%%%%%%%%%%%%%%%%%.@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@..%%%%%%%%%%%%%%%%%%%%.@@@@@@@@.%%%%%%%%%%%%%%%%%%%%.=@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@..%%%%%%%%%%%%%%%%%%%:.@@@@@@.*%%%%%%%%%%%%%%%%%%%..@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@..%%%%%%%%%%%%%%%%%%%..@@@@..%%%%%%%%%%%%%%%%%%%..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..%%%%%%%%%%%%%%%%%%%..@@..%%%%%%%%%%%%%%%%%%%..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..%%%%%%%%%%%%%%%%%%%.:..%%%%%%%%%%%%%%%%%%%..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.%%%%%%%%%%%%%%%%%%%%..%%%%%%%%%%%%%%%%%%%-.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....................@@@@@@@@@@@@@@@@@@@@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..................@@@..@@@@@@@@@@@@@@@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#........................@@@@@@@@@@@@@.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@...............@........+.#@@@@@@@@@.-@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..%%%%%%%%%%%@@+@@....*..%%%@@@%@@..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..%%%%%%%%%%:.@@@@@....%%%%%%+@%..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..%%%%%%%%%%%*@@@@@%%%%%%%%%%%..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@..@..................................................@.:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@..@...@@..@@@..%@...@@@@@..@@@@@...@@@@@@...@@@@@....@.:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@..@...%@..@@@@.=@..@%......@@.....@@....@@.@@........@.:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@..@...%@..@@.@@.@.=@....@@...@@@..@%....@@.@@........@.:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@..@...@@..@@..@@@..@@...@@.....@@.@@-..@@...@@.......@.:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@..@.........................@@.......@=........=.....@.:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@......................%%%%%%%%%%%%.....................:@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.:%%%%%%%%%%..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.%%%%%%%%%*.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.%%%%%%%%.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@.%%%%%%.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@=.%%%%.@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..%%..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@....@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@..@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
                @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            %}
            % Made with https://www.asciiart.eu/image-to-ascii
            % Interpert this as the the gear ratio or its inverse can be an integer to enable a 
            %   step up then step down (or vice versa), as implied by the figure scaling
            if (~is_int(max([N_2, N_3])/min([N_2, N_3])) || ~is_int(max([N_4, N_5])/min([N_4, N_5])))
                continue
            end

            % "Have gear set 2-3 have the larger gear ratio of the two gear sets"
            % if (N_3 / N_2 <= N_5 / N_4)
            if (max([N_2 N_3]) / min([N_2 N_3]) <= max([N_4 N_5]) / min([N_4 N_5]))
                continue;
            end


            OD_2 = get_outer_diam(N_2, diametral_pitch); % [in]
            OD_3 = get_outer_diam(N_3, diametral_pitch); % [in]
            OD_4 = get_outer_diam(N_4, diametral_pitch); % [in]
            OD_5 = get_outer_diam(N_5, diametral_pitch); % [in]

            pitch_dia_2 = get_pitch_diam(N_2, diametral_pitch); % [in]
            pitch_dia_3 = get_pitch_diam(N_3, diametral_pitch); % [in]
            pitch_dia_4 = get_pitch_diam(N_4, diametral_pitch); % [in]
            pitch_dia_5 = get_pitch_diam(N_5, diametral_pitch); % [in]


            % Ensure input, output shafts are aligned
            % Ignore division by two to get pitch radius, as it cancels 
            if (~is_zero((pitch_dia_2 + pitch_dia_3) - (pitch_dia_4 + pitch_dia_5)))
                continue;
            end
    
            shaft_dist = pitch_dia_2/2 + pitch_dia_3/2;

            % It's possible the pitch diameter sums are equivalent, but one set of gears 
            %   has teeth that protrude out further
            working_gearbox_size = max([OD_2 OD_5])/2 + shaft_dist + max([OD_3 OD_4])/2 + ...
                2*(clearance + wall_thickness); % [in]

            % Update gear train and gearbox size if new values are better (or if it wasn't 
            %   initialized yet).
            if (isnan(gearbox_size) || gearbox_size > working_gearbox_size)
                gear_train = [N_2, N_3, N_4, N_5]; % [teeth]
                gearbox_size = working_gearbox_size; % [in]
            end
        end
    end
end

%% Display

% Save values for reference
% MATLAB is rather oddly designed, this just assigns N_2 to the first value of gear_train, N_2 
%   to the 2nd, etc.
[N_2, N_3, N_4, N_5] = feval(@(x) x{:}, num2cell(gear_train));

print_and_validate("N_2", int16(N_2), "teeth", MIN_TEETH_NUM, MAX_TEETH_NUM);
print_and_validate("N_3", int16(N_3), "teeth", MIN_TEETH_NUM, MAX_TEETH_NUM);
print_and_validate("N_4", int16(N_4), "teeth", MIN_TEETH_NUM, MAX_TEETH_NUM);
print_and_validate("N_5", int16(N_5), "teeth", MIN_TEETH_NUM, MAX_TEETH_NUM);
print_and_validate("Y", gearbox_size, "in", 25, 60);


%% Functions

% get_outer_diam
%
% Finds the bounding outside diameter of a gear (including teeth length). A spur gear specified 
%   by the function inputs will fit within a cylidrical envolope with this diameter.
%
% Inputs:
%   - num_teeth -- a scalar double that represents the number of teeth the gear has
%   - diametral_pitch -- a scalar double that represents the number of teeth per length of the gear 
%       - must have the same length units as num_teeth
% 
% Outputs:
%   - OD -- a scalar double that represents the outer diameter of the gear
%
function OD = get_outer_diam(num_teeth, diametral_pitch)
    % See Slide 5 of 34 in 330_S25_Lecture25_Gears
    % [teeth / (teeth/in) ]
    pitch_diameter = num_teeth / diametral_pitch; % [length]

    % See Table 13-1
    addendum = 1/diametral_pitch; % [in], teeth unit magically disappears

    OD = pitch_diameter + addendum*2; % [length]
end

% get_pitch_diam
%
% Finds the bounding outside diameter of a gear (including teeth length). A spur gear specified 
%   by the function inputs will fit within a cylidrical envolope with this diameter.
%
% Inputs:
%   - num_teeth -- a scalar double that represents the number of teeth the gear has
%   - diametral_pitch -- a scalar double that represents the number of teeth per length of the gear 
%       - must have the same length units as num_teeth
% 
% Outputs:
%   - pitch_dia -- a scalar double that represents the outer diameter of the gear
%
function pitch_dia = get_pitch_diam(num_teeth, diametral_pitch)
    % See Slide 5 of 34 in 330_S25_Lecture25_Gears
    % [teeth / (teeth/in) ]
    pitch_dia = num_teeth / diametral_pitch; % [length]
end

% is_int
% 
% Tells you whether a given input is an integer. Requires the variable epsilon to be global.
%
% Inputs:
%   - val -- a double, can be an array
%       - must be less than the maximum integer int16 can store
%
% Outputs:
%   - is_it -- a logical value describing whether val is an integer or not
function is_it = is_int(val)
    global epsilon

    % Ensure everything is cast to double so we don't get any weird conversions back to int and lose 
    %   decimal places.
    % MATLAB doesn't match standard practice of every sensical programming language of truncating 
    %   when casting to integers (it rounds), so need abs as well
    is_it = abs(double(val) - double(int16(val))) < epsilon;
end

% is_zero
% 
% Tells you whether a given input is a zero. Requires the variable epsilon to be global.
%
% Inputs:
%   - val -- a double, can be an array
%       - must be less than the maximum integer int16 can store
%
% Outputs:
%   - is_it -- a logical value describing whether val is a zero
function is_it = is_zero(val)
    global epsilon

    % Ensure everything is cast to double so we don't get any weird conversions back to int and lose 
    %   decimal places.
    % MATLAB doesn't match standard practice of every sensical programming language of truncating 
    %   when casting to integers (it rounds), so need abs as well
    is_it = abs(val) < epsilon;
end


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
%       - Must be 3 characters or less
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
        fprintf("%-3s = %-7s %-5s\n", name, val, units);
        return;
    end
    
    if (contains(class(val), "int"))
        % done instead of %i to avoid leading zero's
        fprintf("%-3s = %-7s %-5s\n", name, string(val), units);
    else
        % format as double
        fprintf("%-3s = %-7.3f %-5s\n", name, val, units);
    end
    
    % Check if range is NaN, if not, proceed with range checking
    if (~sum(isnan([range_low range_high])))
        assert(val >= range_low && val <= range_high, name + " is outside valid range given in"...
            + " problem statement. Manual intervention is required. Good luck.");
    end
end
