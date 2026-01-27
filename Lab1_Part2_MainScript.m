%{
Author(s): Daniel Ghoreyshi, Graeme Appel
Course: ASEN 3802 -- Structures Lab -- Task 2
Goal: Complete Applications Regarding Beam Deflection and force analysis

%}


%% Housekeeping
clc
clear
close all


%% Task 1, Part 1 -- Analysis of Experimental Data

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


%%  Task 2, Part 1 -- Case 1,2,3 Scatter Plots and Best Fits
 
% F0

figure()
scatter(data_Case1_Load, data_Case1_F0, "g", "Linewidth", 0.5);
hold on
scatter(data_Case2_Load, data_Case2_F0, "r", "Linewidth", 0.5);
scatter(data_Case3_Load, data_Case3_F0, "m", "Linewidth", 0.5);

% Find Best fits for each case using polyfit/val

[coeffs_Case1_F0, error_Case1_F0] = polyfit(data_Case1_Load, data_Case1_F0, 1);
R_squared.Case1.F0 = error_Case1_F0.rsquared; % obtain R^2
best_fit_Case1_F0 = polyval(coeffs_Case1_F0, data_Case1_Load);
plot(data_Case1_Load,  best_fit_Case1_F0,"g","Linewidth", 0.5);

[coeffs_Case2_F0, error_Case2_F0] = polyfit(data_Case2_Load, data_Case2_F0, 1);
R_squared.Case2.F0 = error_Case2_F0.rsquared; % obtain R^2
best_fit_Case2_F0 = polyval(coeffs_Case2_F0, data_Case2_Load);
plot(data_Case2_Load,  best_fit_Case2_F0,"r","Linewidth", 0.5);

[coeffs_Case3_F0, error_Case3_F0] = polyfit(data_Case3_Load, data_Case3_F0, 1);
R_squared.Case3.F0 = error_Case3_F0.rsquared; % obtain R^2
best_fit_Case3_F0 = polyval(coeffs_Case3_F0, data_Case3_Load);
plot(data_Case3_Load,  best_fit_Case3_F0,"m","Linewidth", 0.5);

legend("Case 1 Scatter", "Case 2 Scatter", "Case 3 Scatter", ...
    "Case 1 Best Fit","Case 2 Best Fit","Case 3 Best Fit", "Location","best")
title("F0 vs. Applied Loads")
xlabel("Applied load [lbs]")
ylabel("F0 [lbf]")




% F1

figure()
scatter(data_Case1_Load, data_Case1_F1, "g", "Linewidth", 0.5);
hold on
scatter(data_Case2_Load, data_Case2_F1, "r", "Linewidth", 0.5);
scatter(data_Case3_Load, data_Case3_F1, "m", "Linewidth", 0.5);

% Find Best fits for each case using polyfit/val

[coeffs_Case1_F1, error_Case1_F1] = polyfit(data_Case1_Load, data_Case1_F1, 1);
R_squared.Case1.F1 = error_Case1_F1.rsquared; % obtain R^2
best_fit_Case1_F1 = polyval(coeffs_Case1_F1, data_Case1_Load);
plot(data_Case1_Load,  best_fit_Case1_F1,"g","Linewidth", 0.5);

[coeffs_Case2_F1, error_Case2_F1] = polyfit(data_Case2_Load, data_Case2_F1, 1);
R_squared.Case2.F1 = error_Case2_F1.rsquared; % obtain R^2
best_fit_Case2_F1 = polyval(coeffs_Case2_F1, data_Case2_Load);
plot(data_Case2_Load,  best_fit_Case2_F1,"r","Linewidth", 0.5);

[coeffs_Case3_F1, error_Case3_F1] = polyfit(data_Case3_Load, data_Case3_F1, 1);
R_squared.Case3.F1 = error_Case3_F1.rsquared; % obtain R^2
best_fit_Case3_F1 = polyval(coeffs_Case3_F1, data_Case3_Load);
plot(data_Case3_Load,  best_fit_Case3_F1,"m","Linewidth", 0.5);

legend("Case 1 Scatter", "Case 2 Scatter", "Case 3 Scatter", ...
    "Case 1 Best Fit","Case 2 Best Fit","Case 3 Best Fit", "Location","best")
title("F1 vs. Applied Loads")
xlabel("Applied load [lbs]")
ylabel("F1 [lbf]")


% F2

figure()
scatter(data_Case1_Load, data_Case1_F2, "g", "Linewidth", 0.5);
hold on
scatter(data_Case2_Load, data_Case2_F2, "r", "Linewidth", 0.5);
scatter(data_Case3_Load, data_Case3_F2, "m", "Linewidth", 0.5);

% Find Best fits for each case using polyfit/val

[coeffs_Case1_F2, error_Case1_F2] = polyfit(data_Case1_Load, data_Case1_F2, 1);
R_squared.Case1.F2 = error_Case1_F2.rsquared; % obtain R^2
best_fit_Case1_F2 = polyval(coeffs_Case1_F2, data_Case1_Load);
plot(data_Case1_Load,  best_fit_Case1_F2,"g","Linewidth", 0.5);

[coeffs_Case2_F2, error_Case2_F2] = polyfit(data_Case2_Load, data_Case2_F2, 1);
R_squared.Case2.F2 = error_Case2_F2.rsquared; % obtain R^2
best_fit_Case2_F2 = polyval(coeffs_Case2_F2, data_Case2_Load);
plot(data_Case2_Load,  best_fit_Case2_F2,"r","Linewidth", 0.5);

[coeffs_Case3_F2, error_Case3_F2] = polyfit(data_Case3_Load, data_Case3_F2, 1);
R_squared.Case3.F2 = error_Case3_F2.rsquared; % obtain R^2
best_fit_Case3_F2 = polyval(coeffs_Case3_F2, data_Case3_Load);
plot(data_Case3_Load,  best_fit_Case3_F2,"m","Linewidth", 0.5);

legend("Case 1 Scatter", "Case 2 Scatter", "Case 3 Scatter", ...
    "Case 1 Best Fit","Case 2 Best Fit","Case 3 Best Fit", "Location","best")
title("F2 vs. Applied Loads")
xlabel("Applied load [lbs]")
ylabel("F2 [lbf]")



% F3D

figure()
scatter(data_Case1_Load, data_Case1_F3D, "g", "Linewidth", 0.5);
hold on
scatter(data_Case2_Load, data_Case2_F3D, "r", "Linewidth", 0.5);
scatter(data_Case3_Load, data_Case3_F3D, "m", "Linewidth", 0.5);

% Find Best fits for each case using polyfit/val

[coeffs_Case1_F3D, error_Case1_F3D] = polyfit(data_Case1_Load, data_Case1_F3D, 1);
R_squared.Case1.F3D = error_Case1_F3D.rsquared; % obtain R^2
best_fit_Case1_F3D = polyval(coeffs_Case1_F3D, data_Case1_Load);
plot(data_Case1_Load,  best_fit_Case1_F3D,"g","Linewidth", 0.5);

[coeffs_Case2_F3D, error_Case2_F3D] = polyfit(data_Case2_Load, data_Case2_F3D, 1);
R_squared.Case2.F3D = error_Case2_F3D.rsquared; % obtain R^2
best_fit_Case2_F3D = polyval(coeffs_Case2_F3D, data_Case2_Load);
plot(data_Case2_Load,  best_fit_Case2_F3D,"r","Linewidth", 0.5);

[coeffs_Case3_F3D, error_Case3_F3D] = polyfit(data_Case3_Load, data_Case3_F3D, 1);
R_squared.Case3.F3D = error_Case3_F3D.rsquared; % obtain R^2
best_fit_Case3_F3D = polyval(coeffs_Case3_F3D, data_Case3_Load);
plot(data_Case3_Load,  best_fit_Case3_F3D,"m","Linewidth", 0.5);

legend("Case 1 Scatter", "Case 2 Scatter", "Case 3 Scatter", ...
    "Case 1 Best Fit","Case 2 Best Fit","Case 3 Best Fit", "Location","best")
title("F3D vs. Applied Loads")
xlabel("Applied load [lbs]")
ylabel("F3D [lbf]")


% LVDT

figure()
scatter(data_Case1_Load, data_Case1_LVDT, "g", "Linewidth", 0.5);
hold on
scatter(data_Case2_Load, data_Case2_LVDT, "r", "Linewidth", 0.5);
scatter(data_Case3_Load, data_Case3_LVDT, "m", "Linewidth", 0.5);

% Find Best fits for each case using polyfit/val

[coeffs_Case1_LVDT, error_Case1_LVDT] = polyfit(data_Case1_Load, data_Case1_LVDT, 1);
R_squared.Case1.LVDT = error_Case1_LVDT.rsquared; % obtain R^2
best_fit_Case1_LVDT = polyval(coeffs_Case1_LVDT, data_Case1_Load);
plot(data_Case1_Load,  best_fit_Case1_LVDT,"g","Linewidth", 0.5);

[coeffs_Case2_LVDT, error_Case2_LVDT] = polyfit(data_Case2_Load, data_Case2_LVDT, 1);
R_squared.Case2.LVDT = error_Case2_LVDT.rsquared; % obtain R^2
best_fit_Case2_LVDT = polyval(coeffs_Case2_LVDT, data_Case2_Load);
plot(data_Case2_Load,  best_fit_Case2_LVDT,"r","Linewidth", 0.5);

[coeffs_Case3_LVDT, error_Case3_LVDT] = polyfit(data_Case3_Load, data_Case3_LVDT, 1);
R_squared.Case3.LVDT = error_Case3_LVDT.rsquared; % obtain R^2
best_fit_Case3_LVDT = polyval(coeffs_Case3_LVDT, data_Case3_Load);
plot(data_Case3_Load,  best_fit_Case3_LVDT,"m","Linewidth", 0.5);

legend("Case 1 Scatter", "Case 2 Scatter", "Case 3 Scatter", ...
    "Case 1 Best Fit","Case 2 Best Fit","Case 3 Best Fit", "Location","best")
title("LVDT vs. Applied Loads")
xlabel("Applied load [lbs]")
ylabel("LVDT [in]")






