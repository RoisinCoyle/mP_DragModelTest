%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DistanceAssumption
% Date created: 17.05.22
% Date last mostified: 20.05.22
% Purpose: To test the assumption that the distance travelled until the
% particles reach terminal settling vleocity is not significantly different
% than the distance travelled in that time if the particle was always at
% terminal settling velocity
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data
% <<<<<<<<<<<<<<<<
clear

Table_Dio_SA= readtable("./DragModelsTest/Output/20220517/DioguardiOutputVM_SA.txt", "Delimiter", ",");
Table_Dio_Proj= readtable("./DragModelsTest/Output/20220517/DioguardiOutputVM_Proj.txt", "Delimiter", ",");

Table_BB_SA= readtable("./DragModelsTest/Output/20220517/BagheriOutputVM_SA.txt", "Delimiter", ",");
Table_BB_Proj= readtable("./DragModelsTest/Output/20220517/BagheriOutputVM_Proj.txt", "Delimiter", ",");

Table_Zhang_SA= readtable("./DragModelsTest/Output/20220517/ZhangOutputVM_SA.txt", "Delimiter", ",");
Table_Zhang_Proj= readtable("./DragModelsTest/Output/20220517/ZhangOutputVM_Proj.txt", "Delimiter", ",");

Table_Stokes_SA= readtable("./DragModelsTest/Output/20220517/StokesOutputVM_SA.txt", "Delimiter", ",");
Table_Stokes_Proj= readtable("./DragModelsTest/Output/20220517/StokesOutputVM_Proj.txt", "Delimiter", ",");

%% Extract data required from tables
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

DioProj_WCalc = table2array(Table_Dio_Proj(:, "Wt"));
DioProj_time = table2array(Table_Dio_Proj(:, "Time"));
DioProj_Dist = table2array(Table_Dio_Proj(:, "Distance"));

DioSA_WCalc = table2array(Table_Dio_SA(:, "Wt"));
DioSA_time = table2array(Table_Dio_SA(:, "Time"));
DioSA_Dist = table2array(Table_Dio_SA(:, "Distance"));

BBProj_WCalc = table2array(Table_BB_Proj(:, "Wt_Calc"));
BBProj_time = table2array(Table_BB_Proj(:, "Time"));
BBProj_Dist = table2array(Table_BB_Proj(:, "Distance"));

BBSA_WCalc = table2array(Table_BB_SA(:, "Wt_Calc"));
BBSA_time = table2array(Table_BB_SA(:, "Time"));
BBSA_Dist = table2array(Table_BB_SA(:, "Distance"));

ZCProj_WCalc = table2array(Table_Zhang_Proj(:, "Wt"));
ZCProj_time = table2array(Table_Zhang_Proj(:, "Time"));
ZCProj_Dist = table2array(Table_Zhang_Proj(:, "Distance"));

ZCSA_WCalc = table2array(Table_Zhang_SA(:, "Wt"));
ZCSA_time = table2array(Table_Zhang_SA(:, "Time"));
ZCSA_Dist = table2array(Table_Zhang_SA(:, "Distance"));

StokesProj_WCalc = table2array(Table_Stokes_Proj(:, "Wt"));
StokesProj_time = table2array(Table_Stokes_Proj(:, "Time"));
StokesProj_Dist = table2array(Table_Stokes_Proj(:, "Distance"));

StokesSA_WCalc = table2array(Table_Stokes_SA(:, "Wt"));
StokesSA_time = table2array(Table_Stokes_SA(:, "Time"));
StokesSA_Dist = table2array(Table_Stokes_SA(:, "Distance"));

ESD = table2array(Table_Dio_SA(:, "ESD"));  %ESD will be same for all tables

%% Perform calculation
% <<<<<<<<<<<<<<<<<<<<<<

DistConst_Dio_Proj = zeros(140, 1);
DistConst_Dio_SA = zeros(140, 1);
DistConst_BB_Proj = zeros(140, 1);
DistConst_BB_SA = zeros(140, 1);
DistConst_ZC_Proj = zeros(140, 1);
DistConst_ZC_SA = zeros(140, 1);
DistConst_Stokes_Proj = zeros(140, 1);
DistConst_Stokes_SA = zeros(140, 1);

for i=1:140
    DistConst_Dio_Proj(i) = DioProj_WCalc(i) * DioProj_time(i);
    DistConst_Dio_SA(i) = DioSA_WCalc(i) * DioSA_time(i);
    DistConst_BB_Proj(i) = BBProj_WCalc(i) * BBProj_time(i);
    DistConst_BB_SA(i) = BBSA_WCalc(i) * BBSA_time(i);
    DistConst_ZC_Proj(i) = ZCProj_WCalc(i) * ZCProj_time(i);
    DistConst_ZC_SA(i) = ZCSA_WCalc(i) * ZCSA_time(i);
    DistConst_Stokes_Proj(i) = StokesProj_WCalc(i) * StokesProj_time(i);
    DistConst_Stokes_SA(i) = StokesSA_WCalc(i) * StokesSA_time(i);
end

    Diff_DioProj = zeros(140, 1);
    Diff_DioSA = zeros(140, 1);
    Diff_BBProj = zeros(140, 1);
    Diff_BBSA = zeros(140, 1);
    Diff_ZCProj = zeros(140, 1);
    Diff_ZCSA = zeros(140, 1);
    Diff_StokesProj = zeros(140, 1);
    Diff_StokesSA = zeros(140, 1);

for i = 1:140
    Diff_DioProj(i) = DistConst_Dio_Proj(i) - DioProj_Dist(i);
    Diff_DioSA(i) = DistConst_Dio_SA(i) - DioSA_Dist(i)  ;
    Diff_BBProj(i) = DistConst_BB_Proj(i) - BBProj_Dist(i) ;
    Diff_BBSA(i) = DistConst_BB_SA(i) - BBSA_Dist(i);
    Diff_ZCProj(i) = DistConst_ZC_Proj(i) - ZCProj_Dist(i);
    Diff_ZCSA(i) = DistConst_ZC_SA(i) - ZCSA_Dist(i);
    Diff_StokesProj(i) = DistConst_Stokes_Proj(i) - StokesProj_Dist(i);
    Diff_StokesSA(i) = DistConst_Stokes_SA(i) - StokesSA_Dist(i);
end

%% Plot of dist const vs dist calc

%% Dioguardi

% Fit linear model through intercept:
lm_DioProj = fitlm(DioProj_Dist, DistConst_Dio_Proj', 'y~-1+x1');
m_DioProj = lm_DioProj.Coefficients.Estimate(1);
fitY_DioProj = zeros(140, 1);
%plot(lm_DioProj);
% Generate data using linear model:
n1=[max(DistConst_Dio_Proj), max(DioProj_Dist)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
for i=1:140
    fitY_DioProj(i) = m_DioProj * nVal(i);
end

% Plot data
subplot(1, 2, 1)
plot(DioProj_Dist, DistConst_Dio_Proj, 'ob')
xlabel('Modelled distance (m)')
ylabel('Equivalent distance at constant velocity (m)')
title('Dioguardi Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_DioProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx', m_DioProj), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

% Fit linear model through intercept:
lm_DioSA = fitlm(DioSA_Dist, DistConst_Dio_SA, 'y~-1+x1');
m_DioSA = lm_DioSA.Coefficients.Estimate(1);
%plot(lm_DioSA);
fitY_DioSA = zeros(140, 1);
n1=[max(DistConst_Dio_SA), max(DioSA_Dist)];
nMax = max(n1);
nVal=linspace(0, nMax, 140);
for i=1:140
    fitY_DioSA(i) = m_DioSA * nVal(i);
end

% Plot data 
subplot(1, 2, 2)
plot(DioSA_Dist, DistConst_Dio_SA, 'ob')
xlabel('Modelled distance (m)')
ylabel('Equivalent distance at constant velocity (m)')
title('Dioguardi Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_DioSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx', m_DioSA), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/Distance/DioguardiVM_DistVsDist.jpg', 'Resolution', 300)

%% Bagheri
% Fit linear model through intercept:
lm_BBProj = fitlm(BBProj_Dist, DistConst_BB_Proj', 'y~-1+x1');
m_BBProj = lm_BBProj.Coefficients.Estimate(1);
fitY_BBProj = zeros(140, 1);
%plot(lm_BBProj);
% Generate data using linear model:
n1=[max(DistConst_BB_Proj), max(BBProj_Dist)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
for i=1:140
    fitY_BBProj(i) = m_BBProj * nVal(i);
end

% Plot data
subplot(1, 2, 1)
plot(BBProj_Dist, DistConst_BB_Proj, 'ob')
xlabel('Modelled distance (m)')
ylabel('Equivalent distance at constant velocity (m)')
title('Bagheri Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_BBProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx', m_BBProj), 'location', 'best');
hold off

% Fit linear model through intercept:
lm_BBSA = fitlm(BBSA_Dist, DistConst_BB_SA, 'y~-1+x1');
m_BBSA = lm_BBSA.Coefficients.Estimate(1);
%plot(lm_BBSA);
fitY_BBSA = zeros(140, 1);
n1=[max(DistConst_BB_SA), max(BBSA_Dist)];
nMax = max(n1);
nVal=linspace(0, nMax, 140);
for i=1:140
    fitY_BBSA(i) = m_BBSA * nVal(i);
end

% Plot data 
subplot(1, 2, 2)
plot(BBSA_Dist, DistConst_BB_SA, 'ob')
xlabel('Modelled distance (m)')
ylabel('Equivalent distance at constant velocity (m)')
title('Bagheri Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_BBSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx', m_BBSA), 'location', 'best');
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/Distance/BagheriVM_DistVsDist.jpg', 'Resolution', 300)

%% Zhang
% Fit linear model through intercept:
lm_ZCProj = fitlm(ZCProj_Dist, DistConst_ZC_Proj', 'y~-1+x1');
m_ZCProj = lm_ZCProj.Coefficients.Estimate(1);
fitY_ZCProj = zeros(140, 1);
%plot(lm_ZCProj);
% Generate data using linear model:
n1=[max(DistConst_ZC_Proj), max(ZCProj_Dist)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
for i=1:140
    fitY_ZCProj(i) = m_ZCProj * nVal(i);
end

% Plot data
subplot(1, 2, 1)
plot(ZCProj_Dist, DistConst_ZC_Proj, 'ob')
xlabel('Modelled distance (m)')
ylabel('Equivalent distance at constant velocity (m)')
title('Zhang Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_ZCProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx', m_ZCProj), 'location', 'best');
hold off

% Fit linear model through intercept:
lm_ZCSA = fitlm(ZCSA_Dist, DistConst_ZC_SA, 'y~-1+x1');
m_ZCSA = lm_ZCSA.Coefficients.Estimate(1);
%plot(lm_ZCSA);
fitY_ZCSA = zeros(140, 1);
n1=[max(DistConst_ZC_SA), max(ZCSA_Dist)];
nMax = max(n1);
nVal=linspace(0, nMax, 140);
for i=1:140
    fitY_ZCSA(i) = m_ZCSA * nVal(i);
end

% Plot data 
subplot(1, 2, 2)
plot(ZCSA_Dist, DistConst_ZC_SA, 'ob')
xlabel('Modelled distance (m)')
ylabel('Equivalent distance at constant velocity (m)')
title('Zhang Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_ZCSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx', m_ZCSA), 'location', 'best');
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/Distance/ZhangVM_DistVsDist.jpg', 'Resolution', 300)

%% Stokes

% Fit linear model through intercept:
lm_StokesProj = fitlm(StokesProj_Dist, DistConst_Stokes_Proj', 'y~-1+x1');
m_StokesProj = lm_StokesProj.Coefficients.Estimate(1);
fitY_StokesProj = zeros(140, 1);
%plot(lm_StokesProj);
% Generate data using linear model:
n1=[max(DistConst_Stokes_Proj), max(StokesProj_Dist)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 140);
for i=1:140
    fitY_StokesProj(i) = m_StokesProj * nVal(i);
end

% Plot data
subplot(1, 2, 1)
plot(StokesProj_Dist, DistConst_Stokes_Proj, 'ob')
xlabel('Modelled distance (m)')
ylabel('Equivalent distance at constant velocity (m)')
title('Stokes Model, Projection Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_StokesProj, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx', m_StokesProj), 'location', 'best');
hold off

% Fit linear model through intercept:
lm_StokesSA = fitlm(StokesSA_Dist, DistConst_Stokes_SA, 'y~-1+x1');
m_StokesSA = lm_StokesSA.Coefficients.Estimate(1);
%plot(lm_StokesSA);
fitY_StokesSA = zeros(140, 1);
n1=[max(DistConst_Stokes_SA), max(StokesSA_Dist)];
nMax = max(n1);
nVal=linspace(0, nMax, 140);
for i=1:140
    fitY_StokesSA(i) = m_StokesSA * nVal(i);
end

% Plot data 
subplot(1, 2, 2)
plot(StokesSA_Dist, DistConst_Stokes_SA, 'ob')
xlabel('Modelled distance (m)')
ylabel('Equivalent distance at constant velocity (m)')
title('Stokes Model, Surface Area')
hold on
plot(nVal, nVal, '--r')
plot(nVal, fitY_StokesSA, '-g')
legend('Data', 'y=x', sprintf('y=%2.4fx', m_StokesSA), 'location', 'best');
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/Distance/StokesVM_DistVsDist.jpg', 'Resolution', 300)

%% Quick plot to check normality

%% A) Histogram

subplot(4, 2, 1)
histogram(Diff_DioSA, 10)
title("Dioguardi SA")
ylabel("Count")

subplot(4, 2, 2)
histogram(Diff_DioProj, 10)
title("Dioguardi Proj")
ylabel("Count")

subplot(4, 2, 3)
histogram(Diff_BBSA, 10)
title("Bagheri SA")
ylabel("Count")

subplot(4, 2, 4)
histogram(Diff_BBProj)
title("Bagheri Proj")
ylabel("Count")

subplot(4, 2, 5)
histogram(Diff_ZCSA, 10)
title("Zhang SA")
ylabel("Count")

subplot(4, 2, 6)
histogram(Diff_ZCProj, 10)
title("Zhang Proj")
ylabel("Count")

subplot(4, 2, 7)
histogram(Diff_StokesSA, 10)
title("Stokes SA")
ylabel("Count")

subplot(4, 2, 8)
histogram(Diff_StokesProj, 10)
title("Stokes Proj")
ylabel("Count")

%% B) Normal probability plot

subplot(4, 2, 1)
normplot(Diff_DioSA)
title('Normal Probability Plot: Dioguardi SA')

subplot(4, 2, 2)
normplot(Diff_DioProj)
title('', 'Normal Probability Plot: Dioguardi Proj' )

subplot(4, 2, 3)
normplot(Diff_BBSA)
title('Normal Probability Plot: Bagheri SA')

subplot(4, 2, 4)
normplot(Diff_BBProj)
title('Normal Probability Plot: Bagheri Proj')

subplot(4, 2, 5)
normplot(Diff_ZCSA)
title('Normal Probability Plot: Zhang SA')

subplot(4, 2, 6)
normplot(Diff_ZCProj)
title('Normal Probability Plot: Zhang Proj')

subplot(4, 2, 7)
normplot(Diff_StokesSA)
title('Normal Probability Plot: Stokes SA')

subplot(4, 2, 8)
normplot(Diff_StokesProj)
title('Normal Probability Plot: Stokes Proj')


%% Statistical test
% =====================
% Null hypothesis: The average of the difference in distances is equal to
% zero

% Data isn't normally distributed, therefore can't conduct valid t-test.

% Other statistics:

Average = mean(Diff_DioSA);
StandDev = std(Diff_DioSA);
AvA = zeros(140, 1);
MedA = zeros(140, 1);
ModA = zeros (140, 1);
for i=1:140
    AvA(i) = Average;
    MedA(i) = mode(Diff_DioSA);
    ModA(i) = median(Diff_DioSA);
end

plot([1:140], Diff_DioSA, 'ob')
hold on
plot([1:140], AvA, '--r')
plot([1:140], MedA, '--k')
plot([1:140], ModA, '--g')
%% Save results to a new table
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% A) Dio Proj
% ===============
Results_Dio = zeros(140, 2);

for i=1:140
    Results_Dio(i, 1) = DistConst_Dio_Proj(i);
    Results_Dio(i, 2) = Diff_DioProj(i);
end 

Table_Dio = array2table(Results_Dio, "VariableNames", ...
    {'DistanceAtConstantWt', 'DifferenceInDistance'});

Output_Dio_Proj = [Table_Dio_Proj Table_Dio];
Output_Dio_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Output_Dio_Proj, './DragModelsTest/Output/20220517/DioOutput_ProjDistVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Output_Dio_Proj, './DragModelsTest/Output/20220517/DioOutput_ProjDistVM.xls', 'WriteRowNames', true);

%% B) Dio SA
% ===============
Results_Dio = zeros(140, 2);

for i=1:140
    Results_Dio(i, 1) = DistConst_Dio_SA(i);
    Results_Dio(i, 2) = Diff_DioSA(i);
end 

Table_Dio = array2table(Results_Dio, "VariableNames", ...
    {'DistanceAtConstantWt', 'DifferenceInDistance'});

Output_Dio_SA = [Table_Dio_SA Table_Dio];
Output_Dio_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Output_Dio_SA, './DragModelsTest/Output/20220517/DioOutput_DistSAVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Output_Dio_SA, './DragModelsTest/Output/20220517/DioOutput_DistSAVM.xls', 'WriteRowNames', true);

%% C) BB Proj
% ===============
Results_BB = zeros(140, 2);

for i=1:140
    Results_BB(i, 1) = DistConst_BB_Proj(i);
    Results_BB(i, 2) = Diff_BBProj(i);
end 

Table_BB = array2table(Results_BB, "VariableNames", ...
    {'DistanceAtConstantWt', 'DifferenceInDistance'});

Output_BB_Proj = [Table_BB_Proj Table_BB];
Output_BB_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Output_BB_Proj, './DragModelsTest/Output/20220517/BagheriOutput_ProjDistVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Output_BB_Proj, './DragModelsTest/Output/20220517/BagheriOutput_ProjDistVM.xls', 'WriteRowNames', true);

%% D) BB SA
% ===============
Results_BB = zeros(140, 2);

for i=1:140
    Results_BB(i, 1) = DistConst_BB_SA(i);
    Results_BB(i, 2) = Diff_BBSA(i);
end 

Table_BB = array2table(Results_BB, "VariableNames", ...
    {'DistanceAtConstantWt', 'DifferenceInDistance'});

Output_BB_SA = [Table_BB_SA Table_BB];
Output_BB_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Output_BB_SA, './DragModelsTest/Output/20220517/BagheriOutput_DistSAVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Output_BB_SA, './DragModelsTest/Output/20220517/BagheriOutput_DistSAVM.xls', 'WriteRowNames', true);

%% E) Zhang Proj
% ===============
Results_ZC = zeros(140, 2);

for i=1:140
    Results_ZC(i, 1) = DistConst_ZC_Proj(i);
    Results_ZC(i, 2) = Diff_ZCProj(i);
end 

Table_ZC = array2table(Results_ZC, "VariableNames", ...
    {'DistanceAtConstantWt', 'DifferenceInDistance'});

Output_ZC_Proj = [Table_ZC_Proj Table_ZC];
Output_ZC_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Output_ZC_Proj, './DragModelsTest/Output/20220517/ZhangOutput_ProjDistVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Output_ZC_Proj, './DragModelsTest/Output/20220517/ZhangOutput_ProjDistVM.xls', 'WriteRowNames', true);

%% F) Zhang SA
% ===============
Results_ZC = zeros(140, 2);

for i=1:140
    Results_ZC(i, 1) = DistConst_ZC_SA(i);
    Results_ZC(i, 2) = Diff_ZCSA(i);
end 

Table_ZC = array2table(Results_ZC, "VariableNames", ...
    {'DistanceAtConstantWt', 'DifferenceInDistance'});

Output_ZC_SA = [Table_ZC_SA Table_ZC];
Output_ZC_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Output_ZC_SA, './DragModelsTest/Output/20220517/ZhangOutput_DistSAVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Output_ZC_SA, './DragModelsTest/Output/20220517/ZhangOutput_DistSAVM.xls', 'WriteRowNames', true);

%% G) Stokes Proj
% ===============
Results_Stokes = zeros(140, 2);

for i=1:140
    Results_Stokes(i, 1) = DistConst_Stokes_Proj(i);
    Results_Stokes(i, 2) = Diff_StokesProj(i);
end 

Table_Stokes = array2table(Results_Stokes, "VariableNames", ...
    {'DistanceAtConstantWt', 'DifferenceInDistance'});

Output_Stokes_Proj = [Table_Stokes_Proj Table_Stokes];
Output_Stokes_Proj.Properties.VariableNames(1) = {'Shape'};

writetable(Output_Stokes_Proj, './DragModelsTest/Output/20220517/StokesOutput_ProjDistVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Output_Stokes_Proj, './DragModelsTest/Output/20220517/StokesOutput_ProjDistVM.xls', 'WriteRowNames', true);

%% H) Stokes SA
% ===============
Results_Stokes = zeros(140, 2);

for i=1:140
    Results_Stokes(i, 1) = DistConst_Stokes_SA(i);
    Results_Stokes(i, 2) = Diff_StokesSA(i);
end 

Table_Stokes = array2table(Results_Stokes, "VariableNames", ...
    {'DistanceAtConstantWt', 'DifferenceInDistance'});

Output_Stokes_SA = [Table_Stokes_SA Table_Stokes];
Output_Stokes_SA.Properties.VariableNames(1) = {'Shape'};

writetable(Output_Stokes_SA, './DragModelsTest/Output/20220517/StokesOutput_DistSAVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Output_Stokes_SA, './DragModelsTest/Output/20220517/StokesOutput_DistSAVM.xls', 'WriteRowNames', true);
