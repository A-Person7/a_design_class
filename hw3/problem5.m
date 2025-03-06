clear all;

P_i = 4; % [MPa]
r_0 = 200 / 2; % [mm]

% I'm too lazy to solve analytically; use iterative solving
%   Pick a value of x_min
%   Loop through values of r 
%   Determine the maximum shear stress seen for all r
%   Stop once max shear exceeds allowable

tau_allowable = 25; % [MPa]

% % Guess
% x_min = 100;
% % Educated guess
% x_min = 20;
% It's almost uncanny how accurately I know a good guessing value 
% The PHD candidate of guesses
x_min = 17;
decrease_amount = 0.001;


last_valid_x_min = x_min

while true
    x_min = x_min - decrease_amount;

    r_i = r_0 - x_min;

    r = linspace(r_i, r_0, 500);
    sigma_theta = ((r_i^2 * P_i) / (r_0^2 - r_i^2)) * (1 + (r_0 ./ r).^2);
    sigma_r = ((r_i^2 * P_i) / (r_0^2 - r_i^2)) * (1 - (r_0 ./ r).^2);
    % sigma_z = ((r_i^2 * P_i) / (r_0^2 - r_i^2));

    % Since sigma_theta and sigma_r are largest and smallest the principle stresses respectively,
    %   solve for the max shear there. Vectorize to deal with all radius values, just 
    %   in case.
    tau_max = max(0.5*(sigma_theta - sigma_r));

    % Have we found the limit of x_min?
    if (tau_max > tau_allowable)
        fprintf("Minimum x_min is %.3f.\n", last_valid_x_min)
        break
    end

    last_valid_x_min = x_min
end


