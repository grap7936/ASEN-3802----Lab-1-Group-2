%% Housekeeping
clc
clear
close all

%% Task 2, Part 2 -- Determine Location of Off-Centered Load (Case 2)
% Read in Case files for each case
file1 = 'Case File Data/Case 1 - load in center';
data_Case1 = load(file1);

file2 = 'Case File Data/Case 2 data.txt';
data_Case2 = load(file2);

file3 = 'Case File Data/Case 3 data.txt';
data_Case3 = load(file3);


% Extract different columns of data

% Case 1
data_Case1_Load = data_Case1(:,1); % Units [lbs]
data_Case1_F0 = data_Case1(:,2); % Units [lbf]
data_Case1_F1 = data_Case1(:,3); % Units [lbf]
data_Case1_F2 = data_Case1(:,4); % Units [lbf]
data_Case1_F3D = data_Case1(:,5); % Units [lbf]
data_Case1_LVDT = data_Case1(:,6); % Units [in]

% Case 2
data_Case2_Load = data_Case2(:,1); % Units [lbs]
data_Case2_F0 = data_Case2(:,2); % Units [lbf]
data_Case2_F1 = data_Case2(:,3); % Units [lbf]
data_Case2_F2 = data_Case2(:,4); % Units [lbf]
data_Case2_F3D = data_Case2(:,5); % Units [lbf]
data_Case2_LVDT = data_Case2(:,6); % Units [in]

% Case 3
data_Case3_Load = data_Case3(:,1); % Units [lbs]
data_Case3_F0 = data_Case3(:,2); % Units [lbf]
data_Case3_F1 = data_Case3(:,3); % Units [lbf]
data_Case3_F2 = data_Case3(:,4); % Units [lbf]
data_Case3_F3D = data_Case3(:,5); % Units [lbf]
data_Case3_LVDT = data_Case3(:,6); % Units [in]

% Constants (SI Units)
L_total = 16 * 0.250;       % [m]
E = 69e9;                   % [Pa]
I = 2.473 * 10^-6;          % [m^4] 

% 1. Determine location 'a' using Reaction Forces
% Use polyfit to get the slope (RA/P)
RA_Case2 = data_Case2_F0 + data_Case2_F1;
coeff_RA = polyfit(data_Case2_Load, RA_Case2, 1);
slope_RA_P = coeff_RA(1); 

a_measured = L_total * (1 - slope_RA_P); 

% 2. Identify the likely Node
nodes = 0:0.25:L_total;
[~, node_index] = min(abs(nodes - a_measured));
a_predicted = nodes(node_index);

fprintf('Measured load position: %.3f m\n', a_measured);
fprintf('Likely node position: %.3f m (Node %d)\n', a_predicted, node_index-1);

% 3. Compare Mid-span Deflection at P = 50 lbs
P_ref = 50; % [lbs]

% Experimental: use the LVDT slope from polyfit
coeff_LVDT = polyfit(data_Case2_Load, data_Case2_LVDT, 1);
v_exp_50 = coeff_LVDT(1) * P_ref; % Result in [inches]

% Theoretical (using SI, then converting to inches)
P_N = P_ref * 4.44822; % Convert lbs to Newtons
b = L_total - a_predicted;
x = L_total / 2;

% Beam deflection formula for x < a (midspan is at L/2)
v_theory_meters = (P_N * b * x) / (6 * E * I * L_total) * (L_total^2 - b^2 - x^2);
v_theory_50 = v_theory_meters * 39.3701; % Convert m to inches

fprintf('Experimental Mid-span Deflection: %.4f in\n', v_exp_50);
fprintf('Theoretical Mid-span Deflection: %.4f in\n', v_theory_50);