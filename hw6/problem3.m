% Requires Symbol toolbox
%   or a minute of work to solve a simple equation for low cycle fatigue

clear all;

% Both steel cases have same f value as they have the same S_ut value
% 900 MPa ~ 130.53396396 ksi, use slide 12 of 43 in 330_S25_Lecture14_Fatigue1_FullyReversed.pdf
%   as figure title isn't given
f_st = 0.81;

% Case 1

% What do these values mean? What units are they in? Go read the documentation of get_num_cycles
num_1 = get_num_cycles(750, 600, 900, f_st, "1");

% Case 2 
num_2 = get_num_cycles(750, 800, 900, f_st, "2");

% Case 3
num_3 = get_num_cycles(750, 1.1*10^3, 1.2*10^3, 0.85, "3");

% get_num_cycles
%
% Outputs the number of cycles a given fully reversed loading condition can withstand given 
%   an S_e value of 300 MPa.
%
% Checks static failure, low cycle fatigue, and high cycle fatigue. Only works for materials with 
%   an endurance limit of 300 MPa (no aluminum).
%
% This function prints out relevant information as it solves to make it easier to explain and 
%   show work, and make it clear that all cases are checked appropriately.
% 
% This function uses metric units.
%
% Inputs:
%   - sigma_rev
%       - A fully reversed stress experienced, in MPa, in the form of a scalar of type double
%   - S_y
%       - Yield Strength, in MPa, in the form of a scalar of type double
%   - S_ut
%       - Ultimate (Tensile) Strength, in MPa, in the form of a scalar of type double
%   - f
%       - the f value, in the form of a scalar of type double
%   - label
%       - a name to call a given case (Case (label value here) for relevant printouts to keep 
%           things organized.
%       - Can either be of type char (array) (e.g. '1') or a string (e.g. "1")
% Outputs:
%   - num
%       - the number of cycles the part is expected to survive in cases of fatigue, or 1 if 
%           the part experiences static failure
%       - a scalar of type double
% 
% Side Effects:
%   - prints output to the console
function num = get_num_cycles(sigma_rev, S_y, S_ut, f, label) 
    fprintf("\nAnalyzing case %s:\n", label);
    if (sigma_rev >= S_y)
        % Static failure ahh condition
        fprintf("\t- Case %s fails in static failure.\n", label);
        num = 1;

        fprintf("\t- N = %i", round(num, 0));
        fprintf("\n");
        return;
    end
    
    fprintf("\t- Case %s doesn't fail in static fatigue.\n", label);

    % The same for every case
    S_e = 300; % [MPa]
    
    S_f_E3 = f * S_ut; % [MPa]

    % Check for low cycle fatigue
    if (sigma_rev > S_f_E3)
        fprintf("\t- Case %s fails in low cycle fatigue.\n", label);

        % Too lazy to solve by hand, let n be a standin for num to avoid changing variable types
        syms n;
        assume(n, "positive");
    
        eq = sigma_rev == S_ut * n^(log10(f)/3);
        num = double(solve(eq));

        fprintf("\t- N = %i", round(num, 0));
        fprintf("\n");
        return;
    end


    fprintf("\t- Case %s fails in high cycle fatigue.\n", label);

    a = (S_f_E3)^2 / S_e;
    fprintf("\t- a = %5.3f\n", a);
    b = -1/3 * log10(S_f_E3 / S_e);
    fprintf("\t- b = %5.3f\n", b);


    num = (sigma_rev/a)^(1/b);


    fprintf("\t- N = %i", round(num, 0));
    fprintf("\n");
end
