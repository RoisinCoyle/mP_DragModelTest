%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: YuScript: VM
% Date created: 25.04.22
% Date last mostified: 25.04.22
% Purpose: To test the model by Yu satisfies the density and initial
% velocity assumption
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
Density_Dataset = readtable("SettlingVelocity calc\DensityTestTable.txt");

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

%% Yu' method 
% <<<<<<<<<<<<<<<<<
d_dimYu = zeros(54, 1);
wvel_Yu = zeros(54, 1);
CdSph_Yu = zeros(54, 1);
Cd_Yu = zeros(54, 1);

for i=1:54	
    d_dimYu(i) = (((rho_rel(i)*g)/(vis_kin(i)^2.0))^(1.0/3.0))*d_equi(i);
    CdSph_Yu(i) = (432.0/(d_dimYu(i)^3.0))*((1 + 0.022*(d_dimYu(i)^3.0))^0.54)...
                   + (0.47*(1- exp(-0.15*(d_dimYu(i)^0.45))));
    Cd_Yu(i) = CdSph_Yu(i)/(((d_dimYu(i)^-0.25)*(shape_sph(i)^(d_dimYu(i)^0.03))*(CSF(i)^(d_dimYu(i)^0.33)))^0.25);
    wvel_Yu(i) = ((vis_kin(i)*g*rho_rel(i))^(1.0/3.0))*(((4*d_dimYu(i))/(3*Cd_Yu(i)))^(0.5));
end

% Store output in one array
Results_Yu = zeros(54, 4);

for i=1:54
    Results_Yu(i, 1) = d_equi(i);
    Results_Yu(i, 2) = CSF(i);
    Results_Yu(i, 3) = wvel_Yu(i);
    Results_Yu(i, 4) = rho_f(i);
end 

% Removing NaN values
n=0;
for i = 1:54
    if (Results_Yu(i, 3) < 0)
        Results_Yu(i, 3) = 0;
    end
end
%%
Table_Yu = array2table(Results_Yu, "VariableNames", ...
    {'ESD', 'CSF', 'Wt', 'rho_f'});

Table_Yu = [Density_Dataset.ID, Density_Dataset.Shape, array2table(rho_p), Table_Yu];

Table_Yu.Properties.VariableNames(1) = {'ID'};
Table_Yu.Properties.VariableNames(2) = {'Shape'};
Table_Yu.Properties.VariableNames(3) = {'ParticleDensity'};

writetable(Table_Yu, './DragModelsTest/Output/20220517/Density/YuOutputDensity.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Yu, './DragModelsTest/Output/20220517/Density/YuOutputDensity.xls', 'WriteRowNames', true);

%% Plot Yu output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Yu= readtable("./DragModelsTest/Output/20220517/Density/YuOutputDensity.txt", "Delimiter", ",");

%% A1) wt against Density
% =====================

plot(Table_Yu{1:6, 'rho_f'}, Table_Yu{1:6, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Yu{7:12, 'rho_f'}, Table_Yu{7:12, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{13:18, 'rho_f'}, Table_Yu{13:18, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{19:24, 'rho_f'}, Table_Yu{19:24, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{25:30, 'rho_f'}, Table_Yu{25:30, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{31:36, 'rho_f'}, Table_Yu{31:36, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{37:42, 'rho_f'}, Table_Yu{37:42, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_Yu{43:48, 'rho_f'}, Table_Yu{43:48, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_Yu{49:54, 'rho_f'}, Table_Yu{49:54, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Fragment', '', '', 'Fibre', '', '', 'Film', '', '', 'NumColumns', 3, 'location', 'southoutside')
title('Yu Model')
ylabel('Terminal settling velocity (m/s)')
xlabel('Fluid Density (kg/m^{3})')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220517/Density/Yu_Density.jpg', 'Resolution', 300)




