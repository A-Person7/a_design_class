clear all;

%% Setup

safety_factor = 2;


% Geometry

% pitch diameter of the spur gear
pitch_dia = 150; % [mm]

% Length from B to C
center_length = 250; % [mm]
% Length from C to D
end_length = 100; % [mm]

% The diameter of the center length of shaft that we're interested in
syms dia; % [mm]

% Loading Conditions

T_A = 340 * 10^3 * [0, 0, 1]; % [N mm]

% Material Properties
S_y = 420; % [MPa]
S_ut = 560; % [MPa]


%% Part a

% Analyzs tress concentrations at the bearing shoulders or at the gear


%% Part b

% Assume that we want infinite life
