%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: YuScript: VM
% Date created: 25.04.22
% Date last mostified: 25.04.22
% Purpose: To test the model by Yu satisfies the density and initial
% velocity assumption
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file

VM_Dataset = readtable("SettlingVelocity calc\VariedDensitySheet.txt");

rho_p = table2array(VM_Dataset(:, "ParticleDensity"));
rho_f = table2array(VM_Dataset(:, "FluidDensity"));
vis_dyn = table2array(VM_Dataset(:, "DynamicViscosity"));
vis_kin = table2array(VM_Dataset(:, "KinematicVisvosity"));

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


% Set up and calculate additional variables:
SA_mP = zeros(18, 1);
SA_EqSph = zeros(18, 1);
Vol_mP = zeros(18, 1);
Mass_mP = zeros(18, 1);
CSF = zeros(18, 1);
rho_rel = zeros(18, 1);
ProjA_ESD = zeros(18, 1);
g=9.81;

for i=1:18
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
d_dimYu = zeros(18, 1);
wvel_Yu = zeros(18, 1);
CdSph_Yu = zeros(18, 1);
Cd_Yu = zeros(18, 1);

for i=1:18	
    d_dimYu(i) = (((rho_rel(i)*g)/(vis_kin(i)^2.0))^(1.0/3.0))*d_equi(i);
    CdSph_Yu(i) = (432.0/(d_dimYu(i)^3.0))*((1 + 0.022*(d_dimYu(i)^3.0))^0.54)...
                   + (0.47*(1- exp(-0.15*(d_dimYu(i)^0.45))));
    Cd_Yu(i) = CdSph_Yu(i)/(((d_dimYu(i)^-0.25)*(shape_sph(i)^(d_dimYu(i)^0.03))*(CSF(i)^(d_dimYu(i)^0.33)))^0.25);
    wvel_Yu(i) = ((vis_kin(i)*g*rho_rel(i))^(1.0/3.0))*(((4*d_dimYu(i))/(3*Cd_Yu(i)))^(0.5));
end

% Store output in one array
Results_Yu = zeros(18, 4);

for i=1:18
    Results_Yu(i, 1) = d_equi(i);
    Results_Yu(i, 2) = CSF(i);
    Results_Yu(i, 3) = wvel_Yu(i);
    Results_Yu(i, 4) = rho_f(i);
end 

Table_Yu = array2table(Results_Yu, "VariableNames", ...
    {'ESD', 'CSF', 'Wt', 'rho_f'});

Table_Yu = [VM_Dataset.Shape Table_Yu];
Table_Yu.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Yu, './DragModelsTest/Output/YuOutputDensity.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Yu, './DragModelsTest/Output/YuOutputDensity.xls', 'WriteRowNames', true);

%% Plot Yu output
% <<<<<<<<<<<<<<<<<<<
Table_Yu= readtable("./DragModelsTest/Output/YuOutputDensity.txt", "Delimiter", ",");

%% A1) wt against Density
% =====================

plot(Table_Yu{1:6, 'rho_f'}, Table_Yu{1:6, 'Wt'}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Yu{7:12, 'rho_f'}, Table_Yu{7:12, 'Wt'}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{13:18, 'rho_f'}, Table_Yu{13:18, 'Wt'}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Fragment', 'Fibre', 'Film', 'NumColumns', 3, 'location', 'southoutside')
title('Yu Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Fluid Density (kg/m^{3})')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/Yu_Density.jpg', 'Resolution', 300)

