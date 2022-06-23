%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: YuScript: VM
% Date created: 23.04.22
% Date last mostified: 23.06.22
% Purpose: To test the implementation of the Yu drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file

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
Percentage_Error_sq = zeros(140, 1);
RMSE_Sum = 0.0;

for i=1:140
    residual(i) = (wvel_Yu(i) - wvel_meas(i));
    Percentage_Error(i) = abs((residual(i) / wvel_meas(i))*100);
    AE_Sum = AE_Sum + Percentage_Error(i);
    Percentage_Error_sq(i) = ((residual(i)/wvel_meas(i))^2)*100;
    RMSE_Sum = RMSE_Sum + Percentage_Error_sq(i);
end

AE = AE_Sum/140;
RMSE = sqrt(RMSE_Sum/140);

% B) Fragments
residual_F3= zeros(80, 1);
Percentage_Error_F3 = zeros(80, 1);
AE_Sum_F3 = 0.0;
Percentage_Error_sq_F3= zeros(80, 1);
RMSE_Sum_F3 = 0.0;

for i=1:80
    residual_F3(i) = (wvel_Yu(i) - wvel_meas(i));
    Percentage_Error_F3(i) = abs((residual_F3(i) / wvel_meas(i))*100);
    AE_Sum_F3 = AE_Sum_F3 + Percentage_Error_F3(i);
    Percentage_Error_sq_F3(i) = ((residual_F3(i)/wvel_meas(i))^2)*100;
    RMSE_Sum_F3 = RMSE_Sum_F3 + Percentage_Error_sq_F3(i);
end

AE_F3 = AE_Sum_F3/80;
RMSE_F3 = sqrt(RMSE_Sum_F3/80);

% C) Fibres 
residual_F2= zeros(20, 1);
Percentage_Error_F2 = zeros(20, 1);
AE_Sum_F2 = 0.0;
Percentage_Error_sq_F2= zeros(20, 1);
RMSE_Sum_F2 = 0.0;

for i=81:100
    residual_F2(i) = (wvel_Yu(i) - wvel_meas(i));
    Percentage_Error_F2(i) = abs((residual_F2(i) / wvel_meas(i))*100);
    AE_Sum_F2 = AE_Sum_F2 + Percentage_Error_F2(i);
    Percentage_Error_sq_F2(i) = ((residual_F2(i)/wvel_meas(i))^2)*100;
    RMSE_Sum_F2 = RMSE_Sum_F2 + Percentage_Error_sq_F2(i);
end

AE_F2 = AE_Sum_F2/20;
RMSE_F2 = sqrt(RMSE_Sum_F2/20);

% D) Films
residual_F1= zeros(40, 1);
Percentage_Error_F1 = zeros(40, 1);
AE_Sum_F1 = 0.0;
Percentage_Error_sq_F1= zeros(40, 1);
RMSE_Sum_F1 = 0.0;

for i=101:140
    residual_F1(i) = (wvel_Yu(i) - wvel_meas(i));
    Percentage_Error_F1(i) = abs((residual_F1(i) / wvel_meas(i))*100);
    AE_Sum_F1 = AE_Sum_F1 + Percentage_Error_F1(i);
    Percentage_Error_sq_F1(i) = ((residual_F1(i)/wvel_meas(i))^2)*100;
    RMSE_Sum_F1 = RMSE_Sum_F1 + Percentage_Error_sq_F1(i);
end

AE_F1 = AE_Sum_F1/40;
RMSE_F1 = sqrt(RMSE_Sum_F1/40);

Error_table_shape = ["All"; "Fragment"; "Fibre"; "Film"];
Error_table_AE = [AE; AE_F3; AE_F2; AE_F1];
Error_table_RMSE = [RMSE; RMSE_F3; RMSE_F2; RMSE_F1];

Error_table = table(Error_table_shape, Error_table_AE, Error_table_RMSE);

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

%% B1) wt against CSF
% ====================

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

%% B2) wt against CSF
% ====================

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

%% C) wt against wt measured
% ============================
Highest(1) = max(Table_Yu.Wt);
Highest(2) = max(Table_Yu.Wt_Meas);
MaxW = max(Highest);
yx=linspace(0, MaxW, 100);

% Method 1: Plot shapes separately
plot(yx, yx, '-k')
hold on
plot(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Yu Model')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW*1.1] )
set(gca,'XLim', [0, MaxW*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_MeasVsCalc.jpg', 'Resolution', 300);

%% D) wt against wt measured with fitted lines
% ===============================================

plot(Table_Yu.('Wt_Meas'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_Yu.('Wt_Meas'), Table_Yu.('Wt'), 1);
px=[min(Table_Yu.('Wt_Meas')) max(Table_Yu.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(0.8*px(2), 0.95*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Yu.("Wt_Meas")\Table_Yu.("Wt");
mx = m*Table_Yu.("Wt_Meas");
plot(Table_Yu.('Wt_Meas'), mx, '-g');
text(0.7*px(2), 0.8*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Yu Model.')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW])
set(gca, 'Xlim', [0, 1.1*MaxW])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_MeasVsCalc_Eqn.jpg', 'Resolution', 300);

%% D2) wt against wt measured using Matlab fitlm function
% ========================================================
%% D2 A) All shapes together
% Fit linear model through the intercept: 
lm_Yu = fitlm(Table_Yu.Wt_Meas, Table_Yu.Wt, 'y~-1+x1');
m_Yu = lm_Yu.Coefficients.Estimate(1);
fitY_Yu = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Yu.Wt), max(Table_Yu.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_Yu.Rsquared.Ordinary(1);
for i=1:140
    fitY_Yu(i) = m_Yu * nVal(i);
end

plot(Table_Yu.Wt_Meas, Table_Yu.Wt, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Yu Model')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_Yu, '--k')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Yu, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_MeasVsCalc_Fit.jpg', 'Resolution', 300);

%% D2 B) Plot all shapes separately with fitted model

% Fit linear model through the intercept: SA
lm_YuSA = fitlm(Table_Yu.Wt_Meas, Table_Yu.Wt, 'y~-1+x1');
m_YuSA = lm_YuSA.Coefficients.Estimate(1);
fitY_YuSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Yu.Wt), max(Table_Yu.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_YuSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_YuSA(i) = m_YuSA * nVal(i);
end

plot(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Yu Model.')
hold on
plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_YuSA, '--k')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_MeasVsCalc_FitShapes.jpg', 'Resolution', 300);

%% D2 C) Plot Fragments only with fitted model

% Fit linear model through the intercept: SA
lm_YuSAF3 = fitlm(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'y~-1+x1');
m_YuSAF3 = lm_YuSAF3.Coefficients.Estimate(1);
fitY_YuSAF3 = zeros(140, 1);
% Generate data using linear model:
n1_F3=[max(Table_Yu{1:80, "Wt"}), max(Table_Yu{1:80, "Wt_Meas"})] ;
nMax_F3 = max(n1_F3);
nVal_F3=linspace(0, nMax_F3, 140);
r_sq_F3 = lm_YuSAF3.Rsquared.Ordinary(1);
for i=1:140
    fitY_YuSAF3(i) = m_YuSAF3 * nVal_F3(i);
end

plot(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Yu Model.')
hold on
plot(nVal_F3, nVal_F3, '-k')
plot(nVal_F3, fitY_YuSAF3, '--b')
legend('Fragments', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuSAF3, r_sq_F3), 'location', 'best');
set(gca,'YLim', [0, nMax_F3*1.1] )
set(gca,'XLim', [0, nMax_F3*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_MeasVsCalc_FitF3.jpg', 'Resolution', 300);

%% D2 D) Plot fibres separately with fitted model

% Fit linear model through the intercept:
lm_YuSAF2 = fitlm(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'y~-1+x1');
m_YuSAF2 = lm_YuSAF2.Coefficients.Estimate(1);
fitY_YuSAF2 = zeros(140, 1);
% Generate data using linear model:
n1_F2=[max(Table_Yu{81:100, "Wt"}), max(Table_Yu{81:100, "Wt_Meas"})] ;
nMax_F2 = max(n1_F2);
nVal_F2=linspace(0, nMax_F2, 140);
r_sq_F2 = lm_YuSAF2.Rsquared.Ordinary(1);
for i=1:140
    fitY_YuSAF2(i) = m_YuSAF2 * nVal_F2(i);
end

plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Yu Model.')
hold on
plot(nVal_F2, nVal_F2, '-k')
plot(nVal_F2, fitY_YuSAF2, '--r')
legend('Fibres', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuSAF2, r_sq_F2), 'location', 'best');
set(gca,'YLim', [0, nMax_F2*1.1] )
set(gca,'XLim', [0, nMax_F2*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_MeasVsCalc_FitF2.jpg', 'Resolution', 300);

%% D2 E) Plot film separately with fitted model

% Fit linear model through the intercept
lm_YuSAF1 = fitlm(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'y~-1+x1');
m_YuSAF1 = lm_YuSAF1.Coefficients.Estimate(1);
fitY_YuSAF1 = zeros(140, 1);
% Generate data using linear model:
n1_F1=[max(Table_Yu{101:140, "Wt"}), max(Table_Yu{101:140, "Wt_Meas"})] ;
nMax_F1 = max(n1_F1);
nVal_F1=linspace(0, nMax_F1, 140);
r_sq_F1 = lm_YuSAF1.Rsquared.Ordinary(1);
for i=1:140
    fitY_YuSAF1(i) = m_YuSAF1 * nVal_F1(i);
end

plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Yu Model.')
hold on
plot(nVal_F1, nVal_F1, '-k')
plot(nVal_F1, fitY_YuSAF1, '--g')
legend('film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuSAF1, r_sq_F1), 'location', 'best');
set(gca,'YLim', [0, nMax_F1*1.1] )
set(gca,'XLim', [0, nMax_F1*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu/YuVM_MeasVsCalc_FitF1.jpg', 'Resolution', 300);

