%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: OutputAllModels
% Date created: 24.05.22
% Date last mostified: 26.05.22
% Purpose: To plot the output of the drag model tests together
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data
clear

Table_Stokes_SA= readtable("./DragModelsTest/Output/20220621/Stokes/StokesOutputVM_SA.txt", "Delimiter", ',');
Table_Stokes_Proj= readtable("./DragModelsTest/Output/20220621/Stokes/StokesOutputVM_Proj.txt", "Delimiter", ",");

Table_Dio_SA= readtable("./DragModelsTest/Output/20220621/Dioguardi/DioguardiOutputVM_SA.txt", "Delimiter", ",");
Table_Dio_Proj= readtable("./DragModelsTest/Output/20220621/Dioguardi/DioguardiOutputVM_Proj.txt", "Delimiter", ",");

Table_BB_SA= readtable("./DragModelsTest/Output/20220621/Bagheri/BagheriOutputVM_SA.txt", "Delimiter", ",");
Table_BB_Proj= readtable("./DragModelsTest/Output/20220621/Bagheri/BagheriOutputVM_Proj.txt", "Delimiter", ",");

Table_Zhang_SA= readtable("./DragModelsTest/Output/20220621/Zhang/ZhangOutputVM_SA.txt", "Delimiter", ",");
Table_Zhang_Proj= readtable("./DragModelsTest/Output/20220621/Zhang/ZhangOutputVM_Proj.txt", "Delimiter", ",");

Table_Frn= readtable("./DragModelsTest/Output/20220621/Francalanci/FrancalanciOutputVM.txt", "Delimiter", ",");
Table_Dietrich_New = readtable('./DragModelsTest/Output/20220621/Dietrich/DietrichOutputVM_NaN.txt', 'Delimiter', ',');
Table_Dietrich = readtable('./DragModelsTest/Output/20220621/Dietrich/DietrichOutputVM.txt', 'Delimiter', ',');
Table_Yu= readtable("./DragModelsTest/Output/20220621/Yu/YuOutputVM.txt", "Delimiter", ",");

%% Plot implicit models

% SA Outputs
subplot(2, 4, 1)
% Stokes
% Fit linear model through the intercept: SA
lm_StokesSA = fitlm(Table_Stokes_SA.Wt_Meas, Table_Stokes_SA.Wt, 'y~-1+x1');
m_StokesSA = lm_StokesSA.Coefficients.Estimate(1);
fitY_StokesSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Stokes_SA.Wt), max(Table_Stokes_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_StokesSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_StokesSA(i) = m_StokesSA * nVal(i);
end

plot(Table_Stokes_SA.Wt_Meas, Table_Stokes_SA.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Stokes Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_StokesSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_StokesSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 4, 2)
% Fit linear model through the intercept: Projected area
lm_StokesProj = fitlm(Table_Stokes_Proj.Wt_Meas, Table_Stokes_Proj.Wt, 'y~-1+x1');
m_StokesProj = lm_StokesProj.Coefficients.Estimate(1);
fitY_StokesProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Stokes_Proj.Wt), max(Table_Stokes_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_StokesProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_StokesProj(i) = m_StokesProj * nVal(i);
end

plot(Table_Stokes_Proj.Wt_Meas, Table_Stokes_Proj.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Stokes Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_StokesProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_StokesProj, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 4, 3)
% Dioguardi
% Fit linear model through the intercept: SA
lm_DioSA = fitlm(Table_Dio_SA.Wt_Meas, Table_Dio_SA.Wt, 'y~-1+x1');
m_DioSA = lm_DioSA.Coefficients.Estimate(1);
fitY_DioSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Dio_SA.Wt), max(Table_Dio_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_DioSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_DioSA(i) = m_DioSA * nVal(i);
end

plot(Table_Dio_SA.Wt_Meas, Table_Dio_SA.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Dioguardi Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_DioSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 4, 4)
% Fit linear model through the intercept: Projected area
lm_DioProj = fitlm(Table_Dio_Proj.Wt_Meas, Table_Dio_Proj.Wt, 'y~-1+x1');
m_DioProj = lm_DioProj.Coefficients.Estimate(1);
fitY_DioProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Dio_Proj.Wt), max(Table_Dio_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_DioProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_DioProj(i) = m_DioProj * nVal(i);
end

plot(Table_Dio_Proj.Wt_Meas, Table_Dio_Proj.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Dioguardi Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_DioProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioProj, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 4, 5)
% Bagheri
% Fit linear model through the intercept: SA
lm_BBSA = fitlm(Table_BB_SA.Wt_Meas, Table_BB_SA.Wt_Calc, 'y~-1+x1');
m_BBSA = lm_BBSA.Coefficients.Estimate(1);
fitY_BBSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_BB_SA.Wt_Calc), max(Table_BB_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_BBSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_BBSA(i) = m_BBSA * nVal(i);
end

plot(Table_BB_SA.Wt_Meas, Table_BB_SA.Wt_Calc, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Bagheri Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_BBSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 4, 6)
% Fit linear model through the intercept: Projected area
lm_BBProj = fitlm(Table_BB_Proj.Wt_Meas, Table_BB_Proj.Wt_Calc, 'y~-1+x1');
m_BBProj = lm_BBProj.Coefficients.Estimate(1);
fitY_BBProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_BB_Proj.Wt_Calc), max(Table_BB_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_BBProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_BBProj(i) = m_BBProj * nVal(i);
end

plot(Table_BB_Proj.Wt_Meas, Table_BB_Proj.Wt_Calc, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Bagheri Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_BBProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBProj, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 4, 7)
% Zhang
% Fit linear model through the intercept: SA
lm_ZCSA = fitlm(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'y~-1+x1');
m_ZCSA = lm_ZCSA.Coefficients.Estimate(1);
fitY_ZCSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_SA.Wt), max(Table_Zhang_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_ZCSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_ZCSA(i) = m_ZCSA * nVal(i);
end

plot(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Zhang Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_ZCSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSA, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 4, 8)
% Fit linear model through the intercept: Projected area
lm_ZCProj = fitlm(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'y~-1+x1');
m_ZCProj = lm_ZCProj.Coefficients.Estimate(1);
fitY_ZCProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_Proj.Wt), max(Table_Zhang_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_ZCProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_ZCProj(i) = m_ZCProj * nVal(i);
end

plot(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Zhang Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_ZCProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProj, r_sq), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/ImplicitResults_MeasVsCalc_FitAll.jpg', 'Resolution', 300);

%% SA Outputs only
subplot(2, 2, 1)
% Stokes
% Fit linear model through the intercept: SA
lm_StokesSA = fitlm(Table_Stokes_SA.Wt_Meas, Table_Stokes_SA.Wt, 'y~-1+x1');
m_StokesSA = lm_StokesSA.Coefficients.Estimate(1);
fitY_StokesSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Stokes_SA.Wt), max(Table_Stokes_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_StokesSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_StokesSA(i) = m_StokesSA * nVal(i);
end

plot(Table_Stokes_SA.Wt_Meas, Table_Stokes_SA.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Stokes Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_StokesSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_StokesSA, r_sq), 'location', 'bestoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 2, 2)
% Dioguardi
% Fit linear model through the intercept: SA
lm_DioSA = fitlm(Table_Dio_SA.Wt_Meas, Table_Dio_SA.Wt, 'y~-1+x1');
m_DioSA = lm_DioSA.Coefficients.Estimate(1);
fitY_DioSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Dio_SA.Wt), max(Table_Dio_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_DioSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_DioSA(i) = m_DioSA * nVal(i);
end

plot(Table_Dio_SA.Wt_Meas, Table_Dio_SA.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Dioguardi Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_DioSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioSA, r_sq), 'location', 'bestoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 2, 3)
% Bagheri
% Fit linear model through the intercept: SA
lm_BBSA = fitlm(Table_BB_SA.Wt_Meas, Table_BB_SA.Wt_Calc, 'y~-1+x1');
m_BBSA = lm_BBSA.Coefficients.Estimate(1);
fitY_BBSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_BB_SA.Wt_Calc), max(Table_BB_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_BBSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_BBSA(i) = m_BBSA * nVal(i);
end

plot(Table_BB_SA.Wt_Meas, Table_BB_SA.Wt_Calc, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Bagheri Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_BBSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBSA, r_sq), 'location', 'bestoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 2, 4)
% Zhang
% Fit linear model through the intercept: SA
lm_ZCSA = fitlm(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'y~-1+x1');
m_ZCSA = lm_ZCSA.Coefficients.Estimate(1);
fitY_ZCSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_SA.Wt), max(Table_Zhang_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_ZCSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_ZCSA(i) = m_ZCSA * nVal(i);
end

plot(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Zhang Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_ZCSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSA, r_sq), 'location', 'bestoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/ImplicitResults_MeasVsCalc_FitSA.jpg', 'Resolution', 300);

%% Projection area outputs
subplot(2, 2, 1)
% Fit linear model through the intercept: Projected area
lm_StokesProj = fitlm(Table_Stokes_Proj.Wt_Meas, Table_Stokes_Proj.Wt, 'y~-1+x1');
m_StokesProj = lm_StokesProj.Coefficients.Estimate(1);
fitY_StokesProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Stokes_Proj.Wt), max(Table_Stokes_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_StokesProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_StokesProj(i) = m_StokesProj * nVal(i);
end

plot(Table_Stokes_Proj.Wt_Meas, Table_Stokes_Proj.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Stokes Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_StokesProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_StokesProj, r_sq), 'location', 'bestoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 2, 2)
% Fit linear model through the intercept: Projected area
lm_DioProj = fitlm(Table_Dio_Proj.Wt_Meas, Table_Dio_Proj.Wt, 'y~-1+x1');
m_DioProj = lm_DioProj.Coefficients.Estimate(1);
fitY_DioProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Dio_Proj.Wt), max(Table_Dio_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_DioProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_DioProj(i) = m_DioProj * nVal(i);
end

plot(Table_Dio_Proj.Wt_Meas, Table_Dio_Proj.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Dioguardi Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_DioProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioProj, r_sq), 'location', 'bestoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 2, 3)
% Fit linear model through the intercept: Projected area
lm_BBProj = fitlm(Table_BB_Proj.Wt_Meas, Table_BB_Proj.Wt_Calc, 'y~-1+x1');
m_BBProj = lm_BBProj.Coefficients.Estimate(1);
fitY_BBProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_BB_Proj.Wt_Calc), max(Table_BB_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_BBProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_BBProj(i) = m_BBProj * nVal(i);
end

plot(Table_BB_Proj.Wt_Meas, Table_BB_Proj.Wt_Calc, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Bagheri Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_BBProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBProj, r_sq), 'location', 'bestoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(2, 2, 4)
% Fit linear model through the intercept: Projected area
lm_ZCProj = fitlm(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'y~-1+x1');
m_ZCProj = lm_ZCProj.Coefficients.Estimate(1);
fitY_ZCProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_Proj.Wt), max(Table_Zhang_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_ZCProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_ZCProj(i) = m_ZCProj * nVal(i);
end

plot(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Zhang Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_ZCProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProj, r_sq), 'location', 'bestoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/ImplicitResults_MeasVsCalc_FitProj.jpg', 'Resolution', 300);

%% Plot explicit models: Linear models

subplot(1, 3, 1)
% Dietrich
% Fit linear model through the intercept: SA
lm_Dietrich = fitlm(Table_Dietrich_New.Wt_Meas, Table_Dietrich_New.Wt, 'y~-1+x1');
m_Dietrich = lm_Dietrich.Coefficients.Estimate(1);
fitY_Dietrich = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Dietrich_New.Wt), max(Table_Dietrich_New.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_Dietrich.Rsquared.Ordinary(1);
for i=1:140
    fitY_Dietrich(i) = m_Dietrich * nVal(i);
end

plot(Table_Dietrich_New.Wt_Meas, Table_Dietrich_New.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Dietrich Model')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_Dietrich, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Dietrich, r_sq), 'location', 'southoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(1, 3, 2)
% Francalanci
% Fit linear model through the intercept: 
lm_Frn = fitlm(Table_Frn.Wt_Meas, Table_Frn.Wt, 'y~-1+x1');
m_Frn = lm_Frn.Coefficients.Estimate(1);
fitY_Frn = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Frn.Wt), max(Table_Frn.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_Frn.Rsquared.Ordinary(1);
for i=1:140
    fitY_Frn(i) = m_Frn * nVal(i);
end

plot(Table_Frn.Wt_Meas, Table_Frn.Wt, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Francalanci Model')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_Frn, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Frn, r_sq), 'location', 'southoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

subplot(1, 3, 3)
% Yu
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
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Yu Model')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_YuSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuSA, r_sq), 'location', 'southoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/ExplicitResults_MeasVsCalc_Fit.jpg', 'Resolution', 300);

%% Plot explicit models: Shapes on w_meas vs w_calc

% Dietrich
subplot(1, 3, 1)
% Dietrich
% Fit linear model through the intercept: SA
lm_Dietrich = fitlm(Table_Dietrich_New.Wt_Meas, Table_Dietrich_New.Wt, 'y~-1+x1');
m_Dietrich = lm_Dietrich.Coefficients.Estimate(1);
fitY_Dietrich = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Dietrich_New.Wt), max(Table_Dietrich_New.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_Dietrich.Rsquared.Ordinary(1);
for i=1:140
    fitY_Dietrich(i) = m_Dietrich * nVal(i);
end

plot(Table_Dietrich{1:80, "Wt_Meas"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Dietrich{81:100, "Wt_Meas"}, Table_Dietrich{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "Wt_Meas"}, Table_Dietrich{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '--r')
plot(nVal, fitY_Dietrich, '-g')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Dietrich Model')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Dietrich, r_sq), 'NumColumns', 2, 'Location', 'southoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

% Francalanci
subplot(1, 3, 2)
% Fit linear model through origin using fitlm function
lm_Frn = fitlm(Table_Frn.Wt_Meas, Table_Frn.Wt, 'y~-1+x1');
m_Frn = lm_Frn.Coefficients.Estimate(1);
fitY_Frn = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Frn.Wt), max(Table_Frn.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_Frn.Rsquared.Ordinary(1);
for i=1:140
    fitY_Frn(i) = m_Frn * nVal(i);
end

plot(Table_Frn{1:80, "Wt_Meas"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Frn{81:100, "Wt_Meas"}, Table_Frn{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "Wt_Meas"}, Table_Frn{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '--r')
plot(nVal, fitY_Frn, '-g')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Francalanci Model')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Frn, r_sq), 'NumColumns', 2, 'Location', 'southoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

% Yu
subplot(1, 3, 3)
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
hold on
plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '--r')
plot(nVal, fitY_YuSA, '-g')
ylabel(sprintf('Estimated settling \n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('Yu Model')
legend('Fragment', 'Fibre', 'Film', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_YuSA, r_sq), 'NumColumns', 2, 'location', 'southoutside');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
axis square
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/ExplicitResults_MeasVsCalc_FitShape.jpg', 'Resolution', 300);

%% Plot all models: Fitted lines, all shapes

subplot(3, 3, 1) %Yu

% Fit linear model through the intercept
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

plot(Table_Yu{1:80, "Wt_Meas"}, Table_Yu{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \r\n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('A. Yu et al. (2022)')
hold on
plot(Table_Yu{81:100, "Wt_Meas"}, Table_Yu{81:100, "Wt"}, 'sr',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{101:140, "Wt_Meas"}, Table_Yu{101:140, "Wt"}, '^g',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_Yu, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', '', '', '', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Yu, r_sq), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

subplot(3, 3, 2) % Dioguardi Proj

% Fit linear model through the intercept: Projected area
lm_DioProj = fitlm(Table_Dio_Proj.Wt_Meas, Table_Dio_Proj.Wt, 'y~-1+x1');
m_DioProj = lm_DioProj.Coefficients.Estimate(1);
fitY_DioProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Dio_Proj.Wt), max(Table_Dio_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq_Proj = lm_DioProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_DioProj(i) = m_DioProj * nVal(i);
end

plot(Table_Dio_Proj{1:80, "Wt_Meas"}, Table_Dio_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \r\n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title(sprintf('B. Dioguardi et al. (2018): \r\n Projected Area of Volume Equivalent Sphere'))
hold on
plot(Table_Dio_Proj{81:100, "Wt_Meas"}, Table_Dio_Proj{81:100, "Wt"}, 'sr',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dio_Proj{101:140, "Wt_Meas"}, Table_Dio_Proj{101:140, "Wt"}, '^g',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_DioProj, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', '', '', '', sprintf('y=%2.4fx, r^{2}=%1.4f', m_DioProj, r_sq_Proj), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

subplot(3, 3, 3) % Bagheri Proj

% Fit linear model through the intercept: Projected area
lm_BBProj = fitlm(Table_BB_Proj.Wt_Meas, Table_BB_Proj.Wt_Calc, 'y~-1+x1');
m_BBProj = lm_BBProj.Coefficients.Estimate(1);
fitY_BBProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_BB_Proj.Wt_Calc), max(Table_BB_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq_Proj = lm_BBProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_BBProj(i) = m_BBProj * nVal(i);
end

plot(Table_BB_Proj{1:80, "Wt_Meas"}, Table_BB_Proj{1:80, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \r\n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title(sprintf('C. Bagheri and Bonadonna (2016): \r\n Projected Area of Volume Equivalent Sphere'))
hold on
plot(Table_BB_Proj{81:100, "Wt_Meas"}, Table_BB_Proj{81:100, "Wt_Calc"}, 'sr',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_BB_Proj{101:140, "Wt_Meas"}, Table_BB_Proj{101:140, "Wt_Calc"}, '^g',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_BBProj, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', '', '', '', sprintf('y=%2.4fx, r^{2}=%1.4f', m_BBProj, r_sq_Proj), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

subplot(3, 3, 4) % Francalanci

% Fit linear model through the intercept: SA
lm_Frn = fitlm(Table_Frn.Wt_Meas, Table_Frn.Wt, 'y~-1+x1');
m_Frn = lm_Frn.Coefficients.Estimate(1);
fitY_Frn = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Frn.Wt), max(Table_Frn.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_Frn.Rsquared.Ordinary(1);
for i=1:140
    fitY_Frn(i) = m_Frn * nVal(i);
end

plot(Table_Frn{1:80, "Wt_Meas"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \r\n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('D. Francalanci et al. (2021)')
hold on
plot(Table_Frn{81:100, "Wt_Meas"}, Table_Frn{81:100, "Wt"}, 'sr',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "Wt_Meas"}, Table_Frn{101:140, "Wt"}, '^g',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_Frn, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', '', '', '', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Frn, r_sq), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

subplot(3, 3, 5) % Zhang Proj

% Fit linear model through the intercept: Projected area
lm_ZCProj = fitlm(Table_Zhang_Proj.Wt_Meas, Table_Zhang_Proj.Wt, 'y~-1+x1');
m_ZCProj = lm_ZCProj.Coefficients.Estimate(1);
fitY_ZCProj = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_Proj.Wt), max(Table_Zhang_Proj.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq_Proj = lm_ZCProj.Rsquared.Ordinary(1);
for i=1:140
    fitY_ZCProj(i) = m_ZCProj * nVal(i);
end

plot(Table_Zhang_Proj{1:80, "Wt_Meas"}, Table_Zhang_Proj{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \r\n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title(sprintf('E. Zhang and Choi (2021): \r\n Projected Area using Max CSA'))
hold on
plot(Table_Zhang_Proj{81:100, "Wt_Meas"}, Table_Zhang_Proj{81:100, "Wt"}, 'sr',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_Proj{101:140, "Wt_Meas"}, Table_Zhang_Proj{101:140, "Wt"}, '^g',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_ZCProj, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', '', '', '', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCProj, r_sq_Proj), '', '', 'location', 'best');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
hold off

subplot(3, 3, 6) % Zhang SA

% Fit linear model through the intercept: SA
lm_ZCSA = fitlm(Table_Zhang_SA.Wt_Meas, Table_Zhang_SA.Wt, 'y~-1+x1');
m_ZCSA = lm_ZCSA.Coefficients.Estimate(1);
fitY_ZCSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Zhang_SA.Wt), max(Table_Zhang_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq_SA = lm_ZCSA.Rsquared.Ordinary(1);
for i=1:140
    fitY_ZCSA(i) = m_ZCSA * nVal(i);
end

plot(Table_Zhang_SA{1:80, "Wt_Meas"}, Table_Zhang_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \r\n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title(sprintf('F. Zhang and Choi (2021): \r\n Particle Surface Area'))
hold on
plot(Table_Zhang_SA{81:100, "Wt_Meas"}, Table_Zhang_SA{81:100, "Wt"}, 'sr',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Zhang_SA{101:140, "Wt_Meas"}, Table_Zhang_SA{101:140, "Wt"}, '^g',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, fitY_ZCSA, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', '', '', '', sprintf('y=%2.4fx, r^{2}=%1.4f', m_ZCSA, r_sq_SA), '', '', 'location', 'best');
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
hold off

subplot(3, 3, 7) % Dietrich

% Fit linear model through the intercept
lm_Dietrich = fitlm(Table_Dietrich_New.Wt_Meas, Table_Dietrich_New.Wt, 'y~-1+x1');
m_Dietrich = lm_Dietrich.Coefficients.Estimate(1);
fitY_Dietrich = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Dietrich_New.Wt), max(Table_Dietrich_New.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
r_sq = lm_Dietrich.Rsquared.Ordinary(1);
for i=1:140
    fitY_Dietrich(i) = m_Dietrich * nVal(i);
end

plot(Table_Dietrich_New{1:41, "Wt_Meas"}, Table_Dietrich_New{1:41, "Wt"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Dietrich{42, "Wt_Meas"}, Table_Dietrich{42, "Wt"}, 'sr',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel(sprintf('Estimated settling \r\n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('G. Dietrich (1982)')
plot(nVal, nVal, '-k')
plot(nVal, fitY_Dietrich, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('', '', '', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Dietrich, r_sq), '', '', 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

subplot(3, 3, 8) % Stokes Proj

% Fit linear model through the intercept: SA
lm_StokesSA = fitlm(Table_Stokes_SA.Wt_Meas, Table_Stokes_SA.Wt, 'y~-1+x1');
m_StokesSA = lm_StokesSA.Coefficients.Estimate(1);
fitY_StokesSA = zeros(140, 1);
% Generate data using linear model:
n1=[max(Table_Stokes_SA.Wt), max(Table_Stokes_SA.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, 1, 1000);
r_sqSA = lm_StokesSA.Rsquared.Ordinary(1);
for i=1:1000
    fitY_StokesSA(i) = m_StokesSA * nVal(i);
end

plot(Table_Stokes_SA{1:80, "Wt_Meas"}, Table_Stokes_SA{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel(sprintf('Estimated settling \r\n velocity (m/s)'))
xlabel('Measured settling velocity (m/s)')
title('H. Stokes (1851): Particle Surface Area')
hold on
plot(Table_Stokes_SA{81:100, "Wt_Meas"}, Table_Stokes_SA{81:100, "Wt"}, 'sr',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Stokes_SA{101:140, "Wt_Meas"}, Table_Stokes_SA{101:140, "Wt"}, '^g',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(nVal, nVal, '-k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
lgnd = legend('Fragment', 'Fibre', 'Film', 'Calculated Velocity = Measured Velocity', 'Measured velocity +/- 30%', '', 'NumColumns', 2);
lgnd.Position(1) = 0.6916;
lgnd.Position(2) = 0.1126;
lgnd.Position(3) = 0.2134;
lgnd.Position(4) = 0.1914;
title(lgnd, 'Key', 'FontWeight', 'bold');
set(gca,'YLim', [0.0007, 1] )
set(gca,'XLim', [0.0007, 1] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

subplot(3, 3, 9)
plot(nVal, nVal, '-k')
hold on
plot(nVal, fitY_StokesSA, '--k', 'LineWidth', 1)
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
set(gca,'YLim', [0.0005, nMax*1.3] )
set(gca,'XLim', [0.0005, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
set(subplot(3,3,9), 'XTick', [], 'YTick', [])
set(subplot(3,3,9), 'XTickLabel', [], 'YTickLabel', [])
set(gca, 'Position', [0.4108, 0.1126, 0.2134, 0.1914])
axis off
lgnd = legend('', sprintf('y=%2.4fx, r^{2}=%1.4f', m_StokesSA, r_sqSA), '', '');
lgnd.Position(1) = 0.5159;
lgnd.Position(2) = 0.1227;
lgnd.Position(3) = 0.1051;
lgnd.Position(4) = 0.0509;

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/AllFit.jpg', 'Resolution', 300);
