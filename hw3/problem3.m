% All forces in pounds, distances in inches, moments in lbf in, stresses in psi 

clear all;

% Sign
sign_weight = [0, -200, 0]; % lbf
sign_cg = 12 * [10-2/2, -4/2-1/2, 0]; % in
% Take wind force per unit area (acting in - z direction), multiply by PROJECTED area
sign_drag = [0, 0, -2] * (2*12 * 4*12); % lbf
sign_cp = sign_cg;

% Post
post_weight = [0, -300, 0]; % lbf
post_cg = 12 * [10/2, 0, 0]; % in
% Take wind force per unit area (acting in - z direction), multiply by PROJECTED area
post_drag = [0, 0, -2] * (10*12 * 1*12); % lbf
post_cp = post_cg;

% From static equilibrium constraint
R_F = -(sign_weight + sign_drag + post_weight + post_drag);
R_M = -(cross(sign_cg, sign_weight) + cross(sign_cp, sign_drag) + ...
    cross(post_cg, post_weight) + cross(post_cp, post_drag));



% Radius
r = 12/2; % in 
% second area moment of inertia for a circle along its centroid
I = pi/4 * r^4; % in^4
% Make use of additive property of integrals and rotational symmetry of circles
J = 2*I; % in^4

% Break reactive moment into a components (bending moments and torsion
T = R_M(1); % lbf in
M_y = R_M(2); % lbf in
M_z = R_M(3); % lbf in


% Maximum tensile stress location is going to occur in the quadrant with +y, +z
% Find angle at which max stress occurs (90 degrees from where the moment points in the projected 
%   circle)
theta = 90 + atan2d(M_y, M_z); % degrees
assert(theta >= 0 & theta < 90, "??")
% direction of interest
dir_interest = [0, sind(theta), cosd(theta)]; 

r_interest = dir_interest * r; % You guessed it, inches

% % Shear force is y and z component of reactive force perpendicular to the radius
% V = norm(cross(dir_interest, R_F)); % lbf
% A = pi * r^2; % in^2


sigma_x = norm([M_y, M_z]) * r / I; % psi
% Shear due to torsion, don't care about direction
tau_radial = T * r / J; % psi
% No direct shear, as it points in direction of bending

% mohr_3d_fcn(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_xz)
% mohr_3d_fcn(sigma_x, 0, 0, 0, tau_radial, 0)

sigma_y = 0;
sigma_z = 0;
tau_xy = tau_radial;; 
tau_yz = 0;
tau_xz = 0;

mohr_3d_fcn(sigma_x, sigma_y, sigma_z, tau_xy, tau_yz, tau_xz)
