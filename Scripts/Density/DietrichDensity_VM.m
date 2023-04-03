%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DietrichScript: VM
% Date created: 23.04.22
% Date last mostified: 02.03.23
% Purpose: To test the implementation of the Dietrich drag model on a range of
%          particle shapes
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
% Van Mekelebeke (2020) DOI: 10.1021/acs.est.9b07378
% ====================================================
Density_Dataset = readtable("SettlingVelocity calc\DensityTestTableNew.txt");

rho_p = table2array(Density_Dataset(:, "ParticleDensity"));
rho_f = table2array(Density_Dataset(:, "FluidDensity"));
vis_dyn = table2array(Density_Dataset(:, "DynamicViscosity"));
vis_kin = table2array(Density_Dataset(:, "KinematicViscosity"));

d_equi = table2array(Density_Dataset(:, "ParticleSize"));
size_a = table2array(Density_Dataset(:, "a"));
size_b = table2array(Density_Dataset(:, "b"));
size_c = table2array(Density_Dataset(:, "c"));
shape = table2array(Density_Dataset(:, "Shape"));

shape_flt = table2array(Density_Dataset(:, "Flatness"));
shape_eln = table2array(Density_Dataset(:, "elongation"));
shape_del = table2array(Density_Dataset(:, "Dellino"));
shape_sph = table2array(Density_Dataset(:, "Sphericity"));
shape_cir = table2array(Density_Dataset(:, "Circularity"));
Powers = table2array(Density_Dataset(:, "Powers"));


% Set up and calculate additional variables:
SA_mP = zeros(54, 1);
SA_EqSph = zeros(54, 1);
Vol_mP = zeros(54, 1);
Mass_mP = zeros(54, 1);
CSF = zeros(54, 1);
rho_rel = zeros(54, 1);
ProjA_ESD = zeros(54, 1);
g=9.81;

for i=1:54
    SA_EqSph(i) = 4.0*pi()*((d_equi(i)/2.0)^2.0);
    SA_mP(i) = SA_EqSph(i)/shape_sph(i);
    Vol_mP(i) = (4/3)*pi()*((d_equi(i)/2.0)^3.0);
    Mass_mP(i) = rho_p(i)*Vol_mP(i);
    CSF(i) = size_c(i)/(sqrt((size_a(i)*size_b(i))));
    rho_rel(i) = (rho_p(i)-rho_f(i))/rho_f(i);
    ProjA_ESD(i) = pi()*(d_equi(i)^2)*0.25;
end

%% Dietrich' method 1
% <<<<<<<<<<<<<<<<<
% Dietrich's model cannot be used when CSF<0.2
% This is not an iterative procedure, it just calculates the terminal velocity.

% Set up variable arrays
d_dim = zeros(54, 1);
wdim_Dietrich = zeros(54, 1);
wt_Dietrich = zeros(54, 1);
R1 = zeros(54,1);
R2 = zeros(54,1);
R3 = zeros(54,1);
g=9.81;

for i=1:54
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
Results_Dietrich = zeros(54, 3);

for i=1:54
    Results_Dietrich(i, 1) = d_equi(i);
    Results_Dietrich(i, 2) = CSF(i);
    Results_Dietrich(i, 3) = wt_Dietrich(i);
end 

Table_Dietrich = array2table(Results_Dietrich, "VariableNames", ...
    {'ESD', 'CSF', 'Wt_Calc'});

Table_Dietrich = [Density_Dataset.Shape Table_Dietrich];
Table_Dietrich.Properties.VariableNames(1) = {'Shape'};

% writetable(Table_Dietrich, './DragModelsTest/Output/20220621/Density/DietrichOutputVM_Den.txt', 'Delimiter', ',', 'WriteRowNames', true);
% writetable(Table_Dietrich, './DragModelsTest/Output/20220621/Density/DietrichOutputVM_Den.xls', 'WriteRowNames', true);

%% Removing NaN values
n=0;
for i = 1:54
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

writetable(Table_Dietrich_NaN, './DragModelsTest/Output/20220621/Density/DietrichOutputVM_NaNDensity.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dietrich_NaN, './DragModelsTest/Output/20220621/Density/DietrichOutputVM_NaNDensity.xls', 'WriteRowNames', true);

%% Plot Dietrich output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Dietrich= readtable("./DragModelsTest/Output/20220621/Density/DietrichOutputVM_Den.txt", "Delimiter", ",");
Table_Dietrich_New = readtable('./DragModelsTest/Output/20220621/Density/DietrichOutputVM_NaNDensity.txt', 'Delimiter', ',');


%% Plot output after removing NaN Values


%% A!) Density Vs Wt Calc
% =================
plot(Table_Dietrich{1:21, "FluidDensity"}, Table_Dietrich{1:21, "Wt_Calc"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Dietrich{22:42, "FluidDensity"}, Table_Dietrich{22:42, "Wt_Calc"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{43:54, "FluidDensity"}, Table_Dietrich{43:54, "Wt_Calc"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dietrich Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel(sprintf('Fluid density (kg/m^{3})'))
hold off

set(gcf, 'WindowState', 'maximized');
%exportgraphics(gcf, './DragModelsTest/Output/2022051De7/Density/DietrichVM_NaN.jpg', 'Resolution', 300)

%% E3) wt against wt measured using Matlab fitlm function
% ========================================================

% Fit linear model through the intercept: SA
lm_Dietrich = fitlm(Table_Dietrich_New.FluidDensty, Table_Dietrich_New.Wt_Calc, 'y~1');
m_Dietrich = lm_Dietrich.Coefficients.Estimate(1);
fitY_Dietrich = zeros(54, 1);
% Generate data using linear model:
n1=[max(Table_Dietrich_New.FluidDensity), max(Table_Dietrich_New.Wt_Calc)] ;
nMax = max(n1);
nVal=linspace(0, nMax, 54);
r_sq = lm_Dietrich.Rsquared.Ordinary(1);
for i=1:54
    fitY_Dietrich(i) = m_Dietrich * nVal(i);
end

plot(Table_Dietrich_New.FluidDensity, Table_Dietrich_New.Wt_Calc, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
ylabel('Estimated settling velocity (m/s)')
xlabel(sprintf('Fluid Density (kg/m^{3})'))
title('Dietrich Model')
hold on
plot(nVal, fitY_Dietrich, '-g')
legend('Data', sprintf('y=%2.4f, r^{2}=%1.4f', m_Dietrich), 'location', 'best');
set(gca,'YLim', [0, nMax*1.1] )
set(gca,'XLim', [0, nMax*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/DietrichVM_MeasVsCalc_Fit.jpg', 'Resolution', 1200);