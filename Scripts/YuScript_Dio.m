%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: YuScript: VM
% Date created: 23.04.22
% Date last mostified: 23.04.22
% Purpose: To test the implementation of the Yu drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file

%% Read in data file

% Dioguardi (2018) DOI: 10.1002/2017JB014926
% ====================================================
Dio_Dataset = readtable("SettlingVelocity calc\DioguardiSIDataSet.txt");

rho_p = table2array(Dio_Dataset(1:200, "ParticleDensity"));
rho_f = table2array(Dio_Dataset(1:200, "FluidDensity"));
vis_dyn = table2array(Dio_Dataset(1:200, "DynamicViscosity"));

d_equi = table2array(Dio_Dataset(1:200, "ParticleSize"));
size_a = table2array(Dio_Dataset(1:200, "a"));
size_b = table2array(Dio_Dataset(1:200, "b"));
size_c = table2array(Dio_Dataset(1:200, "c"));

shape_flt = table2array(Dio_Dataset(1:200, "flatness"));
shape_eln = table2array(Dio_Dataset(1:200, "elongation"));
shape_del = table2array(Dio_Dataset(1:200, "Dellino"));
shape_sph = table2array(Dio_Dataset(1:200, "Sphericity"));
shape_cir = table2array(Dio_Dataset(1:200, "Circularity"));
Reynolds = table2array(Dio_Dataset(1:200, "Re"));

wvel_meas = table2array(Dio_Dataset(:, "Wmeasured"));

% Set up and calculate additional variables:
SA_mP = zeros(200, 1);
SA_EqSph = zeros(200, 1);
Vol_mP = zeros(200, 1);
Mass_mP = zeros(200, 1);
CSF = zeros(200, 1);
rho_rel = zeros(200, 1);
vis_kin = zeros(200, 1);
ProjA_ESD = zeros(200, 1);
g=9.81;

for i=1:200
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    vis_kin(i) = vis_dyn(i) / rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end


%% Yu' method 
% <<<<<<<<<<<<<<<<<
d_dimYu = zeros(200, 1);
wvel_Yu = zeros(200, 1);
CdSph_Yu = zeros(200, 1);
Cd_Yu = zeros(200, 1);

for i=1:200	
    d_dimYu(i) = (((rho_rel(i)*g)/(vis_kin(i)^2.0))^(1.0/3.0))*d_equi(i);
    CdSph_Yu(i) = (432.0/(d_dimYu(i)^3.0))*((1 + 0.022*(d_dimYu(i)^3.0))^0.54)...
                   + (0.47*(1- exp(-0.15*(d_dimYu(i)^0.45))));
    Cd_Yu(i) = CdSph_Yu(i)/(((d_dimYu(i)^-0.25)*(shape_sph(i)^(d_dimYu(i)^0.03))*(CSF(i)^(d_dimYu(i)^0.33)))^0.25);
    wvel_Yu(i) = ((vis_kin(i)*g*rho_rel(i))^(1.0/3.0))*(((4*d_dimYu(i))/(3*Cd_Yu(i)))^(0.5));
end

% Store output in one array
Results_Yu = zeros(200, 4);

for i=1:200
    Results_Yu(i, 1) = d_equi(i);
    Results_Yu(i, 2) = CSF(i);
    Results_Yu(i, 3) = wvel_Yu(i);
    Results_Yu(i, 4) = wvel_meas(i);
end 

Table_Yu = array2table(Results_Yu, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas'});

writetable(Table_Yu, './DragModelsTest/Output/YuOutputDio.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Yu, './DragModelsTest/Output/YuOutputDio.xls', 'WriteRowNames', true);

%% Plot Yu output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Yu= readtable("./DragModelsTest/Output/YuOutputDio.txt", "Delimiter", ",");

%% A1) wt against ESD
% =====================

plot(Table_Yu.('ESD'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Yu.('ESD'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Yu Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/YuDio_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

% Method 1: Shapes Plotted Separately
plot(Table_Yu.('ESD'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Yu.('ESD'), Table_Yu.('Wt'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'NumColumns', 2, 'location', 'southoutside')
title('Yu Model')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/YuDio_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF
% ====================

% Method 1: Plotting all 
plot(Table_Yu.('CSF'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Yu.('CSF'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Yu Model')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/YuDio_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF
% ====================

% Method 1: Shapes Plotted Separately
plot(Table_Yu.('CSF'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Yu.('CSF'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'NumColumns', 2, 'location', 'southoutside')
title('Yu Model')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/YuDio_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured
% ============================
Highest(1) = max(Table_Yu.Wt);
Highest(2) = max(Table_Yu.Wt_Meas);
MaxW = max(Highest);
yx=linspace(0, MaxW, 100);

% Method 1: Plot shapes separately
plot(yx, yx)
hold on
plot(Table_Yu.("Wt_Meas"), Table_Yu.("Wt"), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
title('Yu Model')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('y=x', '', 'location', 'best')
set(gca,'YLim', [0, MaxW*1.1] )
set(gca,'XLim', [0, MaxW*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/YuDio_MeasVsCalc.jpg', 'Resolution', 300);

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
text(0.75*px(2), 0.93*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Yu.("Wt_Meas")\Table_Yu.("Wt");
mx = m*Table_Yu.("Wt_Meas");
plot(Table_Yu.('Wt_Meas'), mx, '-g');
text(0.71*px(2), 0.65*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Yu Model.')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW])
set(gca, 'Xlim', [0, 1.1*MaxW])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/YuDio_MeasVsCalc_Eqn.jpg', 'Resolution', 300);

%% D2) wt against wt measured using Matlab fitlm function
% ========================================================

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

plot(Table_Yu.Wt_Meas, Table_Yu.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Yu Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_YuSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/YuDio_MeasVsCalc_Fit.jpg', 'Resolution', 300);