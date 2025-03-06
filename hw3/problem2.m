% All forces are in lbf, all lengths are in inches, moments are in lbf*in, etc.
% Coordinate system as described in problem

% As described in the problem statement
P = 1000;
F = 500;
M_1 = 100;
M_2 = 75;


% Duck typed languages may be the only valid usecase for (a bastardized version of )
%   Hungarian notion
% Vectorize forces, moments, add positional vectors to describe location of forces
% All force locations with respect to point A (read as from A to location of force X)

% Forces
% Force itself
vec_P = P*[0, -1, 0];
% Location of the force
vec_r_P = [3 + 0.5, 1 + 1 + 6, 0];
vec_F = F*[-1,0,0];
vec_r_F = [0.5+3,1+1+6,0];

% Pure moments
vec_M_1 = M_1*[0, -1, 0];
vec_M_2 = M_2*[-1, 0, 0];


vec_R_A = - (vec_F + vec_P)
vec_M_A = - (cross(vec_r_F, vec_F) + cross(vec_r_P, vec_P) + vec_M_1 + vec_M_2)
