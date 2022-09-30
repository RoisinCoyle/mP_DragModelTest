%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: FrancalanciScript: VM
% Date created: 23.04.22
% Date last mostified: 22.07.22
% Purpose: To test the implementation of the Francalanci drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
% Van Mekelebeke (2020) DOI: 10.1021/acs.est.9b07378
% ====================================================
VM_Dataset = readtable("SettlingVelocity calc\VanMelkebekeSIDataset.txt");

rho_p = table2array(VM_Dataset(:, "ParticleDensity"));
rho_f = table2array(VM_Dataset(:, "FluidDensity"));
vis_dyn = table2array(VM_Dataset(:, "DynamicViscosity"));
vis_kin = table2array(VM_Dataset(:, "KinematicViscosity"));

d_equi = table2array(VM_Dataset(:, "ParticleSize"));
size_a = table2array(VM_Dataset(:, "a"));
size_b = table2array(VM_Dataset(:, "b"));
size_c = table2array(VM_Dataset(:, "c"));
shape = table2array(VM_Dataset(:, "Shape"));

shape_flt = table2array(VM_Dataset(:, "Flatness"));
shape_eln = table2array(VM_Dataset(:, "elongation"));
shape_del = table2array(VM_Dataset(:, "Dellino"));
shape_sph = table2array(VM_Dataset(:, "Sphericity"));
shape_cir = table2array(VM_Dataset(:, "Circularity"));
Reynolds = table2array(VM_Dataset(:, "Re"));
Powers = table2array(VM_Dataset(:, "Powers"));

wvel_meas = table2array(VM_Dataset(:, "Wmeasured"));

% Set up and calculate additional variables:
SA_mP = zeros(140, 1);
SA_EqSph = zeros(140, 1);
Vol_mP = zeros(140, 1);
Mass_mP = zeros(140, 1);
CSF = zeros(140, 1);
rho_rel = zeros(140, 1);
ProjA_ESD = zeros(140, 1);
g=9.81;

for i=1:140
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end

%% Francalanci' method 
% <<<<<<<<<<<<<<<<<
% This is not an interative procedure, it just calculates the terminal velocity.

Dref_Frn = zeros(140, 1);
Ddim_Frn = zeros(140, 1);
wdim_Frn = zeros(140, 1);
wvel_Frn = zeros(140, 1);

for i=1:140
    Dref_Frn(i) = size_a(i)*(CSF(i)^0.34)*((size_b(i)/size_a(i))^0.5);
    Ddim_Frn(i) = Dref_Frn(i)*(((g*rho_rel(i))/(vis_kin(i)^2.0))^(1.0/3.0));
	
	E = size_a(i)*((((size_a(i)^2.0)+(size_b(i)^2.0)+(size_c(i)^2.0))/3.0)^-0.5);
	C1 = 18.0*(E^-0.38);
	C2 = 0.3708*(CSF(i)^-0.1602);
	n = 0.4942*(CSF(i)^-0.059);
	
	wdim_Frn(i) = (Ddim_Frn(i)^2)/(C1+(0.75*C2*(Ddim_Frn(i)^3))^n);
	wvel_Frn(i) = wdim_Frn(i)*((rho_rel(i)*g*vis_kin(i))^(1.0/3.0));
end

% Store output in one array
Results_Frn = zeros(140, 4);

for i=1:140
    Results_Frn(i, 1) = d_equi(i);
    Results_Frn(i, 2) = CSF(i);
    Results_Frn(i, 3) = wvel_Frn(i);
    Results_Frn(i, 4) = wvel_meas(i);
end 

Table_Frn = array2table(Results_Frn, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas'});

Table_Frn = [VM_Dataset.Shape Table_Frn];
Table_Frn.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Frn, './DragModelsTest/Output/20220621/Francalanci/FrancalanciOutputVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Frn, './DragModelsTest/Output/20220621/Francalanci/FrancalanciOutputVM.xls', 'WriteRowNames', true);

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(140, 1);
Percentage_Error = zeros(140, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wvel_Frn(i) - wvel_meas(i));
    Percentage_Error(i) = (residual(i) / wvel_meas(i))*100;
    AE_Sum = AE_Sum + Percentage_Error(i);
    Abs_AE_Sum = Abs_AE_Sum + abs(Percentage_Error(i));
    Percentage_Error_sq(i) = (Percentage_Error(i))^2;
    RMSE_Sum = RMSE_Sum + Percentage_Error_sq(i);
end

AE = AE_Sum/140;
Abs_AE = Abs_AE_Sum/140;
RMSE = sqrt(RMSE_Sum/140);

% B) Fragments
residual_F3= zeros(80, 1);
Percentage_Error_F3 = zeros(80, 1);
AE_Sum_F3 = 0.0;
Abs_AE_Sum_F3 = 0.0;
Percentage_Error_sq_F3= zeros(80, 1);
RMSE_Sum_F3 = 0.0;

for i=1:80
    residual_F3(i) = (wvel_Frn(i) - wvel_meas(i));
    Percentage_Error_F3(i) = (residual_F3(i) / wvel_meas(i))*100;
    AE_Sum_F3 = AE_Sum_F3 + Percentage_Error_F3(i);
    Abs_AE_Sum_F3 = Abs_AE_Sum_F3 + abs(Percentage_Error_F3(i));
    Percentage_Error_sq_F3(i) = (Percentage_Error_F3(i))^2;
    RMSE_Sum_F3 = RMSE_Sum_F3 + Percentage_Error_sq_F3(i);
end

AE_F3 = AE_Sum_F3/80;
Abs_AE_F3 = Abs_AE_Sum_F3/80;
RMSE_F3 = sqrt(RMSE_Sum_F3/80);

% C) Fibres 
residual_F1= zeros(20, 1);
Percentage_Error_F1 = zeros(20, 1);
AE_Sum_F1 = 0.0;
Abs_AE_Sum_F1 = 0.0;
Percentage_Error_sq_F1= zeros(20, 1);
RMSE_Sum_F1 = 0.0;

for i=81:100
    residual_F1(i-80) = (wvel_Frn(i) - wvel_meas(i));
    Percentage_Error_F1(i-80) = (residual_F1(i-80) / wvel_meas(i))*100;
    AE_Sum_F1 = AE_Sum_F1 + Percentage_Error_F1(i-80);
    Abs_AE_Sum_F1 = Abs_AE_Sum_F1 + abs(Percentage_Error_F1(i-80));
    Percentage_Error_sq_F1(i-80) = (Percentage_Error_F1(i-80))^2;
    RMSE_Sum_F1 = RMSE_Sum_F1 + Percentage_Error_sq_F1(i-80);
end

AE_F1 = AE_Sum_F1/20;
Abs_AE_F1 = Abs_AE_Sum_F1/20;
RMSE_F1 = sqrt(RMSE_Sum_F1/20);

% D) Films
residual_F2= zeros(40, 1);
Percentage_Error_F2 = zeros(40, 1);
AE_Sum_F2 = 0.0;
Abs_AE_Sum_F2 = 0.0;
Percentage_Error_sq_F2= zeros(40, 1);
RMSE_Sum_F2 = 0.0;

for i=101:140
    residual_F2(i-100) = (wvel_Frn(i) - wvel_meas(i));
    Percentage_Error_F2(i-100) = (residual_F2(i-100) / wvel_meas(i))*100;
    AE_Sum_F2 = AE_Sum_F2 + Percentage_Error_F2(i-100);
    Abs_AE_Sum_F2 = Abs_AE_Sum_F2 + abs(Percentage_Error_F2(i-100));
    Percentage_Error_sq_F2(i-100) = (Percentage_Error_F2(i-100))^2;
    RMSE_Sum_F2 = RMSE_Sum_F2 + Percentage_Error_sq_F2(i-100);
end

AE_F2 = AE_Sum_F2/40;
Abs_AE_F2 = Abs_AE_Sum_F2/40;
RMSE_F2 = sqrt(RMSE_Sum_F2/40);

Error_table_shape = ["All"; "Fragment"; "Fibre"; "Film"];
Error_table_AE = [AE; AE_F3; AE_F1; AE_F2];
Error_table_Abs_AE = [Abs_AE; Abs_AE_F3; Abs_AE_F1; Abs_AE_F2];
Error_table_RMSE = [RMSE; RMSE_F3; RMSE_F1; RMSE_F2];

Error_table = table(Error_table_shape, Error_table_AE, Error_table_Abs_AE, Error_table_RMSE);

writetable(Error_table, './DragModelsTest/Output/20220621/Francalanci/FrancalanciErrorTableVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Francalanci/FrancalanciErrorTableVM.xls', 'WriteRowNames', true);

%% Plot Francalanci output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Frn= readtable("./DragModelsTest/Output/20220621/Francalanci/FrancalanciOutputVM.txt", "Delimiter", ",");

%% A1) wt against ESD: All
% =========================

plot(Table_Frn.('ESD'), Table_Frn.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Frn.('ESD'), Table_Frn.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Francalanci Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD: Shapes
% =============================

% Method 1: Shapes Plotted Separately
plot(Table_Frn.('ESD'), Table_Frn.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Frn{1:80, "ESD"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{81:100, "ESD"}, Table_Frn{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "ESD"}, Table_Frn{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Francalanci Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF: All
% =========================

% Method 1: Plotting all 
plot(Table_Frn.('CSF'), Table_Frn.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Frn.('CSF'), Table_Frn.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Francalanci Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF: Shapes
% =============================

% Method 1: Shapes Plotted Separately
plot(Table_Frn.('CSF'), Table_Frn.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Frn{1:80, "CSF"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{81:100, "CSF"}, Table_Frn{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "CSF"}, Table_Frn{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Francalanci Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured using Matlab fitlm function
% ========================================================
%% C A) All shapes together
% Fit linear model through the intercept: 
lm_Frn = fitlm(Table_Frn.Wt_Meas, Table_Frn.Wt, 'y~-1+x1');
m_Frn = lm_Frn.Coefficients.Estimate(1);
fitY_Frn = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Frn.Wt), max(Table_Frn.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq = lm_Frn.Rsquared.Ordinary(1);
for i=1:1000
    fitY_Frn(i) = m_Frn * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Frn.Wt_Meas, Table_Frn.Wt, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Francalanci Model, Surface Area')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_Frn, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Frn, r_sq), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.1] )
set(gca,'XLim', [0.003, nMax*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_MeasVsCalc_Fit.jpg', 'Resolution', 300);

%% C B) Plot all shapes separately with fitted model

% Fit linear model through the intercept: SA
lm_Frn = fitlm(Table_Frn.Wt_Meas, Table_Frn.Wt, 'y~-1+x1');
m_Frn = lm_Frn.Coefficients.Estimate(1);
fitY_Frn = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Frn.Wt), max(Table_Frn.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq = lm_Frn.Rsquared.Ordinary(1);
for i=1:1000
    fitY_Frn(i) = m_Frn * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Frn{1:80, "Wt_Meas"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Francalanci Model.')
hold on
plot(Table_Frn{81:100, "Wt_Meas"}, Table_Frn{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "Wt_Meas"}, Table_Frn{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_Frn, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Frn, r_sq), 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_MeasVsCalc_FitShapes.jpg', 'Resolution', 300);

%% C C) Plot Fragments only with fitted model

% Fit linear model through the intercept: SA
lm_FrnF3 = fitlm(Table_Frn{1:80, "Wt_Meas"}, Table_Frn{1:80, "Wt"}, 'y~-1+x1');
m_FrnF3 = lm_FrnF3.Coefficients.Estimate(1);
fitY_FrnF3 = zeros(1000, 1);
% Generate data using linear model:
n1_F3=[max(Table_Frn{1:80, "Wt"}), max(Table_Frn{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0.0001, nMax_F3, 1000);
r_sq_F3 = lm_FrnF3.Rsquared.Ordinary(1);
for i=1:1000
    fitY_FrnF3(i) = m_FrnF3 * nVal_F3(i);
end

subplot(1, 2, 1)
plot(Table_Frn{1:80, "Wt_Meas"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Francalanci et al (2021).')
hold on
plot(nVal_F3, nVal_F3, '-k')
plot(nVal_F3, fitY_FrnF3, '--b')
plot(nVal_F3, 1.3*nVal_F3, ':k')
plot(nVal_F3, 0.7*nVal_F3, ':k')
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_FrnF3, r_sq_F3), 'location', 'best');
set(gca,'YLim', [0.003, nMax_F3*1.1] )
set(gca,'XLim', [0.003, nMax_F3*1.1] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_MeasVsCalc_FitF3.jpg', 'Resolution', 300);

%% C D) Plot fibres separately with fitted model

% Fit linear model through the intercept:
lm_FrnF2 = fitlm(Table_Frn{81:100, "Wt_Meas"}, Table_Frn{81:100, "Wt"}, 'y~-1+x1');
m_FrnF2 = lm_FrnF2.Coefficients.Estimate(1);
fitY_FrnF2 = zeros(1000, 1);
% Generate data using linear model:
n1_F2=[max(Table_Frn{81:100, "Wt"}), max(Table_Frn{81:100, "Wt_Meas"})] ;
nMax_F2 = max(n1_F2);
nVal_F2=linspace(0.0001, nMax_F2, 1000);
r_sq_F2 = lm_FrnF2.Rsquared.Ordinary(1);
for i=1:1000
    fitY_FrnF2(i) = m_FrnF2 * nVal_F2(i);
end

subplot(1, 2, 1)
plot(Table_Frn{81:100, "Wt_Meas"}, Table_Frn{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Francalanci et al (2021).')
hold on
plot(nVal_F2, nVal_F2, '-k')
plot(nVal_F2, fitY_FrnF2, '--r')
plot(nVal_F2, 1.3*nVal_F2, ':k')
plot(nVal_F2, 0.7*nVal_F2, ':k')
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_FrnF2, r_sq_F2), '', '', 'location', 'best');
set(gca,'YLim', [0.008, nMax_F2*1.1] )
set(gca,'XLim', [0.008, nMax_F2*1.1] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_MeasVsCalc_FitF2.jpg', 'Resolution', 300);

%% C E) Plot film separately with fitted model

% Fit linear model through the intercept
lm_FrnF1 = fitlm(Table_Frn{101:140, "Wt_Meas"}, Table_Frn{101:140, "Wt"}, 'y~-1+x1');
m_FrnF1 = lm_FrnF1.Coefficients.Estimate(1);
fitY_FrnF1 = zeros(1000, 1);
% Generate data using linear model:
n1_F1=[max(Table_Frn{101:140, "Wt"}), max(Table_Frn{101:140, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0.0001, nMax_F1, 1000);
r_sq_F1 = lm_FrnF1.Rsquared.Ordinary(1);
for i=1:1000
    fitY_FrnF1(i) = m_FrnF1 * nVal_F1(i);
end

subplot(1, 2, 1)
plot(Table_Frn{101:140, "Wt_Meas"}, Table_Frn{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Francalanci et al (2021).')
hold on
plot(nVal_F1, nVal_F1, '-k')
plot(nVal_F1, fitY_FrnF1, '--g')
plot(nVal_F1, 1.3*nVal_F1, ':k')
plot(nVal_F1, 0.7*nVal_F1, ':k')
legend('Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_FrnF1, r_sq_F1), 'location', 'best');
set(gca,'YLim', [0.004, nMax_F1*1.1] )
set(gca,'XLim', [0.004, nMax_F1*1.1] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Francalanci/FrancalanciVM_MeasVsCalc_FitF1.jpg', 'Resolution', 300);

%% Combine all m and r_sq values into the error table
Error_table = readtable("./DragModelsTest/Output/20220621/Francalanci/FrancalanciErrorTableVM.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All", "Fragment", "Fibre", "Film"];
Var_types = ["double","double"];

Frn_rsq = [r_sq; r_sq_F3; r_sq_F2; r_sq_F1];
Frn_m = [m_Frn; m_FrnF3; m_FrnF2; m_FrnF1];

Frn_Table = array2table([Frn_m Frn_rsq]);
Frn_Table.Properties.VariableNames = Col_names;
Frn_Table.Properties.RowNames = Row_names;

Error_table = [Error_table Frn_Table];

writetable(Error_table, './DragModelsTest/Output/20220621/Francalanci/FrancalanciFinalTableVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Francalanci/FrancalanciFinalTableVM.xls', 'WriteRowNames', true);