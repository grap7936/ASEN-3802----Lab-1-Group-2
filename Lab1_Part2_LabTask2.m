%% ASEN 3802 Lab 1 (Part 2) — Task 2 ONLY (Case 2)
clc; clear; close all;

%% Load Case 2
D2 = load('Case File Data/Case 2 data.txt');

P_lbf   = D2(:,1);
F0_lbf  = D2(:,2);
F1_lbf  = D2(:,3);
F2_lbf  = D2(:,4);
F3D_lbf = D2(:,5);
LVDT_in = D2(:,6);

%% Constants + conversions
L = 16*0.250;     % [m]
E = 69e9;         % [Pa]
I = 2.473e-6;     % [m^4]

lbf2N = 4.4482216152605;
m2in  = 39.370078740157;

% Member properties for Fi prediction
A_in2 = 0.0614;        % [in^2]
d_in  = 4.921;         % [in]
A = A_in2*(0.0254^2);  % [m^2]
c = (d_in/2)*0.0254;   % [m]

%% 1) Experimental load location from reactions
RA_lbf = F0_lbf + F1_lbf;      % left support (per your assumption)
RB_lbf = F2_lbf;               % right support

% Use slope-through-origin for RA/P (more robust than full polyfit for ratios)
sRA = (P_lbf' * RA_lbf) / (P_lbf' * P_lbf);    % ~ RA/P
a_meas_m = L*(1 - sRA);                         % [m] (from left)

% Snap to nearest node
nodes_m = 0:0.25:L;
[~, node_idx] = min(abs(nodes_m - a_meas_m));
a_node_m = nodes_m(node_idx);

%% 2) Experimental midspan deflection at P = 50 lbf (from LVDT fit)
P_ref_lbf = 50;

pV = polyfit(P_lbf, LVDT_in, 1);      % v = mP + b
v_exp_50_in = polyval(pV, P_ref_lbf); % [in]

%% 3) Theoretical midspan deflection at P = 50 lbf (equivalent beam)
x = L/2;                 % midspan
a = a_node_m;            % load location from left [m]
b = L - a;

P_ref_N = P_ref_lbf*lbf2N;

if x <= a
    v_mid_m = (P_ref_N*b*x)/(6*E*I*L) * (L^2 - b^2 - x^2);
else
    v_mid_m = (P_ref_N*a*(L-x))/(6*E*I*L) * (L^2 - a^2 - (L-x)^2);
end

v_th_50_in = v_mid_m*m2in;

pct_err_v = abs(v_exp_50_in - v_th_50_in)/abs(v_th_50_in)*100;

%% 4) Expected internal force Fi vs measured F3D (use slopes, preload-free)
P_N = P_lbf*lbf2N;

% Midspan moment for a single point load
M_mid_Nm = 0.5 .* P_N .* min(a, b);     % [N*m]

% Expected internal force
Fi_expected_lbf = ((A*c/I) .* M_mid_Nm) / lbf2N;

% Preload-correct F3D using intercept (or just compare slopes)
pF = polyfit(P_lbf, F3D_lbf, 1);        % F3D = kP + b
F3D_corrected_lbf = F3D_lbf - pF(2);

% Compare slopes (single-number comparison)
pFi = polyfit(P_lbf, Fi_expected_lbf, 1);
pFc = polyfit(P_lbf, F3D_corrected_lbf, 1);

%% OUTPUTS (Task 2 values only)
fprintf('\n=== Task 2 (Case 2) Results ===\n');
fprintf('a_measured = %.3f m (%.2f in)\n', a_meas_m, a_meas_m*m2in);
fprintf('a_node     = %.3f m (%.2f in)  (Node %d)\n', a_node_m, a_node_m*m2in, node_idx-1);

fprintf('\nMidspan deflection at P = 50 lbf:\n');
fprintf('v_exp(50)  = %.4f in\n', v_exp_50_in);
fprintf('v_th(50)   = %.4f in\n', v_th_50_in);
fprintf('%% error    = %.1f %%\n', pct_err_v);

fprintf('\nInternal force comparison (slope vs load):\n');
fprintf('slope(Fi_expected / P) = %.4f [lbf/lbf]\n', pFi(1));
fprintf('slope(F3D_corrected/P) = %.4f [lbf/lbf]\n', pFc(1));
fprintf('ratio (meas/model)     = %.4f\n', pFc(1)/pFi(1));

fprintf('\nF3D preload estimate (intercept): %.3f lbf\n', pF(2));

% Small summary table at the reference load (50 lbf)
Fi_exp_50 = polyval(pFi, P_ref_lbf);
F3Dcorr_50 = polyval(pFc, P_ref_lbf);


pct_err_F3D = abs(Fi_exp_50 - F3Dcorr_50)/abs(F3Dcorr_50)*100;

fprintf('\nAt P = 50 lbf:\n');
fprintf('Fi_expected(50) ~ %.3f lbf\n', Fi_exp_50);
fprintf('F3D_corrected(50) ~ %.3f lbf\n', F3Dcorr_50);
fprintf('%% error    = %.1f %%\n', pct_err_F3D);