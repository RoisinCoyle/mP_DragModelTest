%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DietrichScript: VM
% Date created: 23.04.22
% Date last mostified: 22.07.22
% Purpose: To test the implementation of the Dietrich drag model on a range of
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

for i=1:140
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end

%% Dietrich
% <<<<<<<<<<<<<<<<<
% Dietrich's model cannot be used when CSF<0.2
% This is not an interative procedure, it just calculates the terminal velocity.

% Set up variable arrays
d_dim = zeros(140, 1);
wdim_Dietrich = zeros(140, 1);
wt_Dietrich = zeros(140, 1);
R1 = zeros(140,1);
R2 = zeros(140,1);
R3 = zeros(140,1);
g=9.81;

for i=1:140
        d_dim(i) = ((rho_p(i) - rho_f(i))*g*(d_equi(i)^3.0))/(rho_f(i)*(vis_kin(i)^2.0));
	
	    R1(i) = -3.76715 + 1.92944*(log10(d_dim(i))) - 0.09815*((log10(d_dim(i)))^2.0) ...
		    -0.00575*((log10(d_dim(i)))^3.0) + 0.00056*((log10(d_dim(i)))^4.0);
	    
        if(CSF(i)>0.15)
            R2(i) = (log10(1-((1-CSF(i)))/0.85)) - ((1-CSF(i))^2.3)*tanh(log10(d_dim(i))-4.6)...
		    + 0.3*(0.5-CSF(i))*((1-CSF(i))^2.0)*(log10(d_dim(i))-4.6);
        else
            R2(i)=nan;
        end
	    R3(i) = (0.65-((CSF(i)/2.83)*(tanh(log10(d_dim(i))-4.6))))^(1+((3.5-Powers(i))/2.5));
	    
	    wdim_Dietrich(i) = R3(i) * (10^(R1(i)+R2(i)));
	    
	    wt_Dietrich(i) = ((wdim_Dietrich(i)*(rho_p(i) - rho_f(i))*g*vis_kin(i))/rho_f(i))^(1.0/3.0);
end

% Store output in one array
Results_Dietrich = zeros(140, 4);

for i=1:140
    Results_Dietrich(i, 1) = d_equi(i);
    Results_Dietrich(i, 2) = CSF(i);
    Results_Dietrich(i, 3) = wt_Dietrich(i);
    Results_Dietrich(i, 4) = wvel_meas(i);
end 

Table_Dietrich = array2table(Results_Dietrich, "VariableNames", ...
    {'ESD', 'CSF', 'Wt', 'Wt_Meas'});

Table_Dietrich = [VM_Dataset.Shape Table_Dietrich];
Table_Dietrich.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Dietrich, './DragModelsTest/Output/20220621/Dietrich/DietrichOutputVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dietrich, './DragModelsTest/Output/20220621/Dietrich/DietrichOutputVM.xls', 'WriteRowNames', true);

%% Removing NaN values
n=0;
for i = 1:140
    if (isnan(Results_Dietrich(i, 3)))
        n=n;
    else
        n = n+1;
        Particle_Num(n) = i;
        Results_NaN_Dietrich(n, :) = Results_Dietrich(i, :);
    end
end

NewShape = Table_Dietrich.Shape(Particle_Num);
NewShape_T = array2table(NewShape);
Table_Dietrich_NaN = array2table(Results_NaN_Dietrich, "VariableNames", ...
    {'ESD', 'CSF', 'Wt', 'Wt_Meas'});

Table_Dietrich_NaN = [NewShape_T Table_Dietrich_NaN];
Table_Dietrich_NaN.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Dietrich_NaN, './DragModelsTest/Output/20220621/Dietrich/DietrichOutputVM_NaN.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dietrich_NaN, './DragModelsTest/Output/20220621/Dietrich/DietrichOutputVM_NaN.xls', 'WriteRowNames', true);

%% Calculate average error and RMSE

% A) All shapes
residual = zeros(42, 1);
Percentage_Error = zeros(42, 1);
AE_Sum = 0.0;
Abs_AE_Sum = 0.0;
Percentage_Error_sq = zeros(42, 1);
RMSE_Sum = 0.0;

for i=1:42
    residual(i) = (Table_Dietrich_NaN.Wt(i)- Table_Dietrich_NaN.Wt_Meas(i));
    Percentage_Error(i) = (residual(i) / Table_Dietrich_NaN.Wt_Meas(i))*100;
    AE_Sum = AE_Sum + Percentage_Error(i);
    Abs_AE_Sum = Abs_AE_Sum + abs(Percentage_Error(i));
    Percentage_Error_sq(i) = (Percentage_Error(i))^2;
    RMSE_Sum = RMSE_Sum + Percentage_Error_sq(i);
end

AE = AE_Sum/42;
Abs_AE = Abs_AE_Sum/42;
RMSE = sqrt(RMSE_Sum/42);

Error_table_shape = ["All"];
Error_table_AE = [AE];
Error_table_Abs_AE = [Abs_AE];
Error_table_RMSE = [RMSE];

Error_table = table(Error_table_shape, Error_table_AE, Error_table_Abs_AE, Error_table_RMSE);

writetable(Error_table, './DragModelsTest/Output/20220621/Dietrich/DietrichErrorTableVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Dietrich/DietrichErrorTableVM.xls', 'WriteRowNames', true);

%% Plot Dietrich output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Dietrich= readtable("./DragModelsTest/Output/20220621/Dietrich/DietrichOutputVM.txt", "Delimiter", ",");
Table_Dietrich_New = readtable('./DragModelsTest/Output/20220621/Dietrich/DietrichOutputVM_NaN.txt', 'Delimiter', ',');

%% A1) wt against ESD
% =====================

plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dietrich Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dietrich/DietrichVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

% Method 1: Shapes Plotted Separately
plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dietrich{1:80, "ESD"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich{81:100, "ESD"}, Table_Dietrich{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "ESD"}, Table_Dietrich{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dietrich Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dietrich/DietrichVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF: All
% =========================

% Method 1: Plotting all 
plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dietrich Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dietrich/DietrichVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF: Shapes
% ============================

% Method 1: Shapes Plotted Separately
plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dietrich{1:80, "CSF"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich{81:100, "CSF"}, Table_Dietrich{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "CSF"}, Table_Dietrich{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dietrich (1982)')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dietrich/DietrichVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) Plot output after removing NaN Values
%% C1) CSF Vs Wt Meas
% =================
subplot(1, 2, 1)
plot(Table_Dietrich_New.('CSF'), Table_Dietrich_New.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dietrich_New{1:41, "CSF"}, Table_Dietrich_New{1:41, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich_New{42, "CSF"}, Table_Dietrich_New{42, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dietrich (1982)')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

subplot(1, 2, 2)
plot(Table_Dietrich_New.('ESD'), Table_Dietrich_New.('Wt_Meas'), 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dietrich_New{1:41, "ESD"}, Table_Dietrich_New{1:41, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich_New{42, "ESD"}, Table_Dietrich_New{42, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dietrich (1982)')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dietrich/DietrichVM_ESDV_CSF_NaN.jpg', 'Resolution', 300)

%% C2) wt against wt measured using Matlab fitlm function alongside CSF
% ======================================================================

% Fit linear model through the intercept
lm_Dietrich = fitlm(Table_Dietrich_New.Wt_Meas, Table_Dietrich_New.Wt, 'y~-1+x1');
m_Dietrich = lm_Dietrich.Coefficients.Estimate(1);
fitY_Dietrich = zeros(1000, 1);
% Generate data using linear model:
n1=[max(Table_Dietrich_New.Wt), max(Table_Dietrich_New.Wt_Meas)] ;
nMax = max(n1);
nVal=linspace(0.0001, nMax*1.3, 1000);
r_sq = lm_Dietrich.Rsquared.Ordinary(1);
for i=1:1000
    fitY_Dietrich(i) = m_Dietrich * nVal(i);
end

subplot(1, 2, 1)
plot(Table_Dietrich_New{1:40, "Wt_Meas"}, Table_Dietrich_New{1:40, "Wt"}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Dietrich{41, "Wt_Meas"}, Table_Dietrich{41, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ylabel('Estimated settling velocity (m/s)')
xlabel('Measured settling velocity (m/s)')
title(sprintf('Dietrich (1982): Estimated terminal settling velocity \n\r against measured terminal settling velocity'))
hold on
plot(nVal, nVal, '-k')
plot(nVal, fitY_Dietrich, '--k')
plot(nVal, 1.3*nVal, ':k')
plot(nVal, 0.7*nVal, ':k')
legend('Fragment', 'Fibre', 'y=x', sprintf('y=%2.4fx, r^{2}=%1.4f', m_Dietrich, r_sq), 'location', 'best');
set(gca,'YLim', [0.003, nMax*1.3] )
set(gca,'XLim', [0.003, nMax*1.3] )
set(gca, 'YScale', 'log')
set(gca, 'XScale', 'log')
hold off

subplot(1, 2, 2)
plot(Table_Dietrich{1:80, 'CSF'}, Table_Dietrich{1:80,'Wt_Meas'}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','b', 'MarkerFaceColor', '[.7, .7, .7]')
hold on
plot(Table_Dietrich{81:100, 'CSF'}, Table_Dietrich{81:100,'Wt_Meas'}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','r', 'MarkerFaceColor', '[.7, .7, .7]')
plot(Table_Dietrich{101:140, 'CSF'}, Table_Dietrich{101:140,'Wt_Meas'}, 'o', ...
    'MarkerSize',5,'MarkerEdgeColor','g', 'MarkerFaceColor', '[.7, .7, .7]')
plot(Table_Dietrich_New{1:41, "CSF"}, Table_Dietrich_New{1:41, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich_New{42, "CSF"}, Table_Dietrich_New{42, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
ymax = max(max(Table_Dietrich_New.Wt_Meas), max(Table_Dietrich_New.Wt));
ydot = linspace(0, ymax*1.1, 100)
xdot = linspace(0.15, 0.15, 100)
plot(xdot, ydot, 'k:')
legend('Measured Wt, Fragment', 'Measured Wt, Fibre', 'Measured Wt, Film', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'CSF=0.15', 'NumColumns', 1, 'location', 'best')
title(sprintf('Dietrich (1982): Corey Shape Factor \n Against Terminal Settling Velocity'))
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
set(gca, 'YLim', [0, ymax*1.1])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Dietrich/DietrichVM_MeasVsCalc_FitCSF.jpg', 'Resolution', 300);

%% Combine all m and r_sq values into the error table
Error_table = readtable("./DragModelsTest/Output/20220621/Dietrich/DietrichErrorTableVM.txt", 'Delimiter', ',', ReadVariableNames=true, ReadRowNames=true);

Col_names = ["m", "r_sq"];
Row_names = ["All"];
Var_types = ["double","double"];

Dietrich_rsq = [r_sq];
Dietrich_m = [m_Dietrich];

Dietrich_Table = array2table([Dietrich_m Dietrich_rsq]);
Dietrich_Table.Properties.VariableNames = Col_names;
Dietrich_Table.Properties.RowNames = Row_names;

Error_table = [Error_table Dietrich_Table];

writetable(Error_table, './DragModelsTest/Output/20220621/Dietrich/DietrichFinalTableVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Error_table, './DragModelsTest/Output/20220621/Dietrich/DietrichFinalTableVM.xls', 'WriteRowNames', true);

