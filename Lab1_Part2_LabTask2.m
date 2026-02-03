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


% Constants (Ensure SI Units)
L = 16 * 0.250;             % Total length (16 bays @ 0.25m) [m]
E = 69e9;                   % Elastic Modulus [Pa]
I = 2.473 *10^4;  

% 1. Determine location 'a' using Reaction Forces
coeff_RA = polyfit(data_Case2_Load, (data_Case2_F0 + data_Case2_F1), 1);
slope_RA = coeff_RA(1); % This is RA/P

a_measured = L * (1 - slope_RA); 

% 2. Identify the likely Node
% Since loads are only at joints, a must be a multiple of 0.25m
nodes = 0:0.25:L;
[~, node_index] = min(abs(nodes - a_measured));
a_predicted = nodes(node_index);

fprintf('Measured load position: %.3f m\n', a_measured);
fprintf('Likely node position: %.3f m (Node %d)\n', a_predicted, node_index-1);

% 3. Compare Mid-span Deflection
% Calculate theoretical v(L/2) using a_predicted and compare to LVDT slope
coeff_LVDT = polyfit(data_Case2_Load, data_Case2_LVDT, 1);
slope_LVDT_exp = coeff_LVDT(1); % [in/lb] -> convert to [m/N] for comparison
