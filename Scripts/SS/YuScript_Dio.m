%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: YuScript: VM
% Date created: 23.04.22
% Date last mostified: 21.07.22
% Purpose: To test the implementation of the Yu drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
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

writetable(Table_Yu, './DragModelsTest/Output/20220621/Yu_Dio/YuOutputDio.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Yu, './DragModelsTest/Output/20220621/Yu_Dio/YuOutputDio.xls', 'WriteRowNames', true);

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(200, 1);
Percentage_Error = zeros(200, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(200, 1);
RMSE_Sum = 0.0;

for i=1:200
    residual(i) = (wvel_Yu(i) - wvel_meas(i));
    Percentage_Error(i) = (residual(i) / wvel_meas(i))*100;
    AE_Sum = AE_Sum + Percentage_Error(i);
    Abs_AE_Sum = Abs_AE_Sum + abs(Percentage_Error(i));
    Percentage_Error_sq(i) = (Percentage_Error(i))^2;
    RMSE_Sum = RMSE_Sum + Percentage_Error_sq(i);
end

AE = AE_Sum/200;
Abs_AE = Abs_AE_Sum/200;
RMSE = sqrt(RMSE_Sum/200);


Error_table_shape = ["All"];
Error_table_AE = [AE];
Error_table_RMSE = [RMSE];

Error_table = table(Error_table_shape, Error_table_AE, Error_table_RMSE);

writetable(Error_table, './DragModelsTest/Output/20220621/Yu_Dio/YuDioErrorTable.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Yu_Dio/YuDioErrorTable.xls', 'WriteRowNames', true);

%% Plot Yu output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Yu= readtable("./DragModelsTest/Output/20220621/Yu_Dio/YuOutputDio.txt", "Delimiter", ",");

%% A1) wt against ESD
% =====================

% Shapes together
subplot(1, 2, 2)
plot(Table_Yu.('ESD'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Yu.('ESD'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Yu (2022): Using Dioguardi Dataset')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu_Dio/YuDio_ESDVsW.jpg', 'Resolution', 300)


%% B1) wt against CSF
% ====================

% Method 1: Plotting all 
subplot(1, 2, 2)
plot(Table_Yu.('CSF'), Table_Yu.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Yu.('CSF'), Table_Yu.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Yu (2022): Using Dioguardi Dataset')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu_Dio/YuDio_CSFVsW.jpg', 'Resolution', 300);

%% D) wt against wt measured using Matlab fitlm function
% ========================================================
%% D A) All shapes together
% Fit linear model through the intercept: 
lm_Yu = fitlm(Table_Yu.Wt_Meas, Table_Yu.Wt, 'y~-1+x1');
m_Yu = lm_Yu.Coefficients.Estimate(1);
fitY_Yu = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Yu.Wt), max(Table_Yu.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.001, 1.1*nMax, 1000);
r_sq = lm_Yu.Rsquared.Ordinary(1);
for i=1:1000
    fitY_Yu(i) = m_Yu * nVal(i);
end

subplot(1, 2, 2)
plot(Table_Yu.Wt_Meas, Table_Yu.Wt, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title('Yu (2022): Using Dioguardi Dataset')
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_Yu, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Yu, r_sq), '+/- 30%', 'location', 'best');
set(gca,'YLim', [0.0015, nMax*1.1] )
set(gca,'XLim', [0.0015, nMax*1.1] )
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Yu_Dio/YuDio_MeasVsCalc_Fit.jpg', 'Resolution', 300);

%% Combine all m and r_sq values into the error table
Error_table = readtable("./DragModelsTest/Output/20220621/Yu_Dio/YuDioErrorTable.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All"];
Var_types = ["double","double"];

Yu_rsq = [r_sq];
Yu_m = [m_Yu];

Yu_Table = array2table([Yu_m Yu_rsq]);
Yu_Table.Properties.VariableNames = Col_names;
Yu_Table.Properties.RowNames = Row_names;

Error_table = [Error_table Yu_Table];

writetable(Error_table, './DragModelsTest/Output/20220621/Yu_Dio/YuDioFinalTable.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Yu_Dio/YuDioFinalTable.xls', 'WriteRowNames', true);