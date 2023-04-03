%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: YuScript: VM
% Date created: 23.04.22
% Date last mostified: 02.03.23
% Purpose: To test the implementation of the Yu drag model on a range of
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

%% Yu' method 
% <<<<<<<<<<<<<<<<<
d_dimYu = zeros(140, 1);
wvel_Yu = zeros(140, 1);
CdSph_Yu = zeros(140, 1);
Cd_Yu = zeros(140, 1);

for i=1:140	
    d_dimYu(i) = (((rho_rel(i)*g)/(vis_kin(i)^2.0))^(1.0/3.0))*d_equi(i);
    CdSph_Yu(i) = (432.0/(d_dimYu(i)^3.0))*((1 + 0.022*(d_dimYu(i)^3.0))^0.54)...
                   + (0.47*(1- exp(-0.15*(d_dimYu(i)^0.45))));
    Cd_Yu(i) = CdSph_Yu(i)/(((d_dimYu(i)^-0.25)*(shape_sph(i)^(d_dimYu(i)^0.03))*(CSF(i)^(d_dimYu(i)^0.33)))^0.25);
    wvel_Yu(i) = ((vis_kin(i)*g*rho_rel(i))^(1.0/3.0))*(((4*d_dimYu(i))/(3*Cd_Yu(i)))^(0.5));
end

% Store output in one array
Results_Yu = zeros(140, 4);

for i=1:140
    Results_Yu(i, 1) = d_equi(i);
    Results_Yu(i, 2) = CSF(i);
    Results_Yu(i, 3) = wvel_Yu(i);
    Results_Yu(i, 4) = wvel_meas(i);
end 

Table_Yu = array2table(Results_Yu, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas'});

Table_Yu = [VM_Dataset.Shape Table_Yu];
Table_Yu.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Yu, './DragModelsTest/Output/20220621/Yu/YuOutputVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Yu, './DragModelsTest/Output/20220621/Yu/YuOutputVM.xls', 'WriteRowNames', true);

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(140, 1);
Percentage_Error = zeros(140, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wvel_Yu(i) - wvel_meas(i));
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
    residual_F3(i) = (wvel_Yu(i) - wvel_meas(i));
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
    residual_F1(i) = (wvel_Yu(i) - wvel_meas(i));
    Percentage_Error_F1(i) = (residual_F1(i) / wvel_meas(i))*100;
    AE_Sum_F1 = AE_Sum_F1 + Percentage_Error_F1(i);
    Abs_AE_Sum_F1 = Abs_AE_Sum_F1 + abs(Percentage_Error_F1(i));
    Percentage_Error_sq_F1(i) = (Percentage_Error_F1(i))^2;
    RMSE_Sum_F1 = RMSE_Sum_F1 + Percentage_Error_sq_F1(i);
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
    residual_F2(i) = (wvel_Yu(i) - wvel_meas(i));
    Percentage_Error_F2(i) = (residual_F2(i) / wvel_meas(i))*100;
    AE_Sum_F2 = AE_Sum_F2 + Percentage_Error_F2(i);
    Abs_AE_Sum_F2 = Abs_AE_Sum_F2 + abs(Percentage_Error_F2(i));
    Percentage_Error_sq_F2(i) = (Percentage_Error_F2(i))^2;
    RMSE_Sum_F2 = RMSE_Sum_F2 + Percentage_Error_sq_F2(i);
end

AE_F2 = AE_Sum_F2/40;
Abs_AE_F2 = Abs_AE_Sum_F2/40;
RMSE_F2 = sqrt(RMSE_Sum_F2/40);

Error_table_shape = ["All"; "Fragment"; "Fibre"; "Film"];
Error_table_AE = [AE; AE_F3; AE_F1; AE_F2];
Error_table_Abs_AE = [Abs_AE; Abs_AE_F3; Abs_AE_F1; Abs_AE_F2];
Error_table_RMSE = [RMSE; RMSE_F3; RMSE_F1; RMSE_F2];

Error_table = table(Error_table_shape, Error_table_AE, Error_table_Abs_AE, Error_table_RMSE);

writetable(Error_table, './DragModelsTest/Output/20220621/Yu/YuErrorTableVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Yu/YuErrorTableVM.xls', 'WriteRowNames', true);

%% Plot Yu output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Yu= readtable("./DragModelsTest/Output/20220621/Yu/YuOutputVM.txt", "Delimiter", ",");

%% A1) wt against ESD
% =====================

plot(Table_Yu.('ESD'), Table_Yu.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Yu.('ESD'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Yu Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

% Method 1: Shapes Plotted Separately
plot(Table_Yu.('ESD'), Table_Yu.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Yu{1:80, "ESD"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{81:100, "ESD"}, Table_Yu{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "ESD"}, Table_Yu{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Yu Model')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF:All
% ========================

% Method 1: Plotting all 
plot(Table_Yu.('CSF'), Table_Yu.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Yu.('CSF'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Yu Model')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF: Shapes
% ============================

% Method 1: Shapes Plotted Separately
plot(Table_Yu.('CSF'), Table_Yu.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Yu{1:80, "CSF"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{81:100, "CSF"}, Table_Yu{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "CSF"}, Table_Yu{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Yu Model')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured using Matlab fitlm function
% ========================================================
%% C A) All shapes together
% Fit linear model through the intercept: 
lm_Yu = fitlm(Table_Yu.Wt_Meas, Table_Yu.Wt, 'y~-1+x1');
m_Yu = lm_Yu.Coefficients.Estimate(1);
fitY_Yu = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Yu.Wt), max(Table_Yu.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq = lm_Yu.Rsquared.Ordinary(1);
for i=1:1000
    fitY_Yu(i) = m_Yu * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Yu.Wt_Meas, Table_Yu.Wt, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP terminal settling velocity to \n\r mP terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Yu et al (2022).'))

hold on
plot(nVal, nVal, '-k', 'LineWidth', 1)
plot(nVal, fitY_Yu, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal, 0.7*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Yu, r_sq), 'y = x +/- 30%', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.1] )
set(gca,'XLim', [0.003, nMax*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Yu/YuVM_MeasVsCalc_Fit.jpg', 'Resolution', 1200);

%% C B) Plot all shapes separately with fitted model

% Fit linear model through the intercept: SA
lm_Yu = fitlm(Table_Yu.Wt_Meas, Table_Yu.Wt, 'y~-1+x1');
m_Yu = lm_Yu.Coefficients.Estimate(1);
fitY_Yu = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Yu.Wt), max(Table_Yu.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax, 1000);
r_sq = lm_Yu.Rsquared.Ordinary(1);
for i=1:1000
    fitY_Yu(i) = m_Yu * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP terminal settling velocity to \n\r mP terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Yu et al (2022).'))
hold on
plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k', 'LineWidth', 1)
plot(nVal, fitY_Yu, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal, 0.7*nVal, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Yu, r_sq), 'y = x +/- 30%', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Yu/YuVM_MeasVsCalc_FitShapes.jpg', 'Resolution', 1200);

%% C C) Plot Fragments only with fitted model

% Fit linear model through the intercept: SA
lm_YuF3 = fitlm(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'y~-1+x1');
m_YuF3 = lm_YuF3.Coefficients.Estimate(1);
fitY_YuF3 = zeros(1000, 1);
% Generate data using linear model:
n1_F3=[max(Table_Yu{1:80, "Wt"}), max(Table_Yu{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0.0001, nMax_F3, 1000);
r_sq_F3 = lm_YuF3.Rsquared.Ordinary(1);
for i=1:1000
    fitY_YuF3(i) = m_YuF3 * nVal_F3(i);
end

subplot(1, 2, 1)
plot(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP fragment terminal settling velocity to \n\r mP fragment terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Yu et al (2022).'))
hold on
plot(nVal_F3, nVal_F3, '-k', 'LineWidth', 1)
plot(nVal_F3, fitY_YuF3, '--b', 'LineWidth', 1)
plot(nVal_F3, 1.3*nVal_F3, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F3, 0.7*nVal_F3, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuF3, r_sq_F3), 'y = x +/- 30%', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F3*1.1] )
set(gca,'XLim', [0.003, nMax_F3*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Yu/YuVM_MeasVsCalc_FitF3.jpg', 'Resolution', 1200);

%% C D) Plot fibres separately with fitted model

% Fit linear model through the intercept:
lm_YuF1 = fitlm(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'y~-1+x1');
m_YuF1 = lm_YuF1.Coefficients.Estimate(1);
fitY_YuF1 = zeros(1000, 1);
% Generate data using linear model:
n1_F1=[max(Table_Yu{81:100, "Wt"}), max(Table_Yu{81:100, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0.0001, nMax_F1, 1000);
r_sq_F1 = lm_YuF1.Rsquared.Ordinary(1);
for i=1:1000
    fitY_YuF1(i) = m_YuF1 * nVal_F1(i);
end

subplot(1, 2, 1)
plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP fibre terminal settling velocity to \n\r mP fibre terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Yu et al (2022).'))
hold on
plot(nVal_F1, nVal_F1, '-k', 'LineWidth', 1)
plot(nVal_F1, fitY_YuF1, '--r', 'LineWidth', 1)
plot(nVal_F1, 1.3*nVal_F1, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F1, 0.7*nVal_F1, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuF1, r_sq_F1), 'y = x +/-30%', 'location', 'best');
set(gca,'YLim', [0.007, nMax_F1*1.1] )
set(gca,'XLim', [0.007, nMax_F1*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Yu/YuVM_MeasVsCalc_FitF1.jpg', 'Resolution', 1200);

%% C E) Plot film separately with fitted model

% Fit linear model through the intercept
lm_YuF2 = fitlm(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'y~-1+x1');
m_YuF2 = lm_YuF2.Coefficients.Estimate(1);
fitY_YuF2 = zeros(1000, 1);
% Generate data using linear model:
n1_F2=[max(Table_Yu{101:140, "Wt"}), max(Table_Yu{101:140, "Wt_Meas"})] ;
nMax_F2 = max(n1_F1);
nVal_F2=linspace(0.0001, nMax_F2, 1000);
r_sq_F2 = lm_YuF2.Rsquared.Ordinary(1);
for i=1:1000
    fitY_YuF2(i) = m_YuF2 * nVal_F2(i);
end

subplot(1, 2, 1)
plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Modelled terminal settling velocity (m/s)')
xlabel('Measured terminal settling velocity (m/s)')
title(sprintf('Graph comparing modelled mP film terminal settling velocity to \n\r mP film terminal settling velocity measured by Van Melkebeke et al (2020).'))
subtitle(sprintf('Model applied: Yu et al (2022).'))
hold on
plot(nVal_F2, nVal_F2, '-k', 'LineWidth', 1)
plot(nVal_F2, fitY_YuF2, '--g', 'LineWidth', 1)
plot(nVal_F2, 1.3*nVal_F2, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
plot(nVal_F2, 0.7*nVal_F2, ':k', 'LineWidth', 1.5, 'Color', [.7 .7 .7])
legend('film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuF2, r_sq_F2), 'y = x +/- 30%', 'location', 'best');
set(gca,'YLim', [0.003, nMax_F1*1.1] )
set(gca,'XLim', [0.003, nMax_F1*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Yu/YuVM_MeasVsCalc_FitF2.jpg', 'Resolution', 1200);

%% Combine all m and r_sq values into the error table
Error_table = readtable("./DragModelsTest/Output/20220621/Yu/YuErrorTableVM.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All", "Fragment", "Fibre", "Film"];
Var_types = ["double","double"];

Yu_rsq = [r_sq; r_sq_F3; r_sq_F1; r_sq_F2];
Yu_m = [m_Yu; m_YuF3; m_YuF1; m_YuF2];

Yu_Table = array2table([Yu_m Yu_rsq]);
Yu_Table.Properties.VariableNames = Col_names;
Yu_Table.Properties.RowNames = Row_names;

Error_table = [Error_table Yu_Table];

writetable(Error_table, './DragModelsTest/Output/20220621/Yu/YuFinalTableVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Yu/YuFinalTableVM.xls', 'WriteRowNames', true);

