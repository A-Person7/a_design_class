clear all;

%% Geometry
w = 1.5; % [in]
l = 6.5; % [in]
r = 0.25; % [in]
d = w - 2*r; % [in]
t = 0.1; % [in]

%% Loading Conditions, Material Properties
M = 40; % [lbf in]
S_ut = 70; % [ksi]



%% Stress

% See Figure A-15-4
c = d/2; % [in]
% c = 0.75
I = t * d^3 / 12; % [in^4]

sigma_0 = M*c/I; % [psi]

fprintf("sigma_0 = %s psi\n\n", addComma(sigma_0));

fprintf("r/d = %3.2f\n", r/d);
fprintf("w/d = %3.2f\n", w/d);

% See Figure A-15-4, will need to be manually updated if r, w, or d change
K_t = 1.6;

fprintf("\nK_t = %.2f\n", K_t);


% Can directly plug in S_ut as it's in [ksi]
sqrt_a = 0.246 - 3.08*10^(-3)*S_ut + 1.51*10^(-5)*S_ut^2 - 2.67*10^(-8)*S_ut^3;
q = 1 / (1 + sqrt_a / sqrt(r));

fprintf("q   = %.2f\n", q);
% Use Fig. 6-26 in George and Waisman 1969 (on slide 28/43 of 
%   330_S25_Lecture14_Fatigue1_FullyReversed.pdf) to verify this q value (expect q ~ 0.8)

K_f = 1 + q*(K_t - 1); % [psi]

fprintf("K_f = %.2f\n", K_f);

sigma_rev = K_f * sigma_0; % [psi]

fprintf("\nsigma_rev = %s psi\n", addComma(sigma_rev));

% Use Java formatting code to get nicer looking numbers (adds commas so 1000 (double) becomes 
%   '1,000' (char)
% MATLAB is kind of secretly Java, some C++, and probably a little left over FORTRAN in a coat
%   anyways
% https://www.mathworks.com/matlabcentral/answers/96131-is-there-a-format-in-matlab-to-display-numbers-such-that-commas-are-automatically-inserted-into-the
function numOut = addComma(numIn)
   import java.text.*
   jf=java.text.DecimalFormat; % comma for thousands, three decimal places
   numOut= char(jf.format(numIn)); % omit "char" if you want a string out
end
