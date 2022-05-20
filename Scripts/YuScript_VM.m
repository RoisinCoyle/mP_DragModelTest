%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: YuScript: VM
% Date created: 23.04.22
% Date last mostified: 23.04.22
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

writetable(Table_Yu, './DragModelsTest/Output/YuOutputVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Yu, './DragModelsTest/Output/YuOutputVM.xls', 'WriteRowNames', true);

%% Plot Yu output
% <<<<<<<<<<<<<<<<<<<
Table_Yu= readtable("./DragModelsTest/Output/YuOutputVM.txt", "Delimiter", ",");

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
exportgraphics(gcf, './DragModelsTest/Output/YuVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

% Method 1: Shapes Plotted Separately
plot(Table_Yu.('ESD'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
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
exportgraphics(gcf, './DragModelsTest/Output/YuVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

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
exportgraphics(gcf, './DragModelsTest/Output/YuVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF
% ====================

% Method 1: Shapes Plotted Separately
plot(Table_Yu.('CSF'), Table_Yu.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
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
exportgraphics(gcf, './DragModelsTest/Output/YuVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured
% ============================
Highest(1) = max(Table_Yu.Wt);
Highest(2) = max(Table_Yu.Wt_Meas);
MaxW = max(Highest);
yx=linspace(0, MaxW, 100);

% Method 1: Plot shapes separately
plot(yx, yx)
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
exportgraphics(gcf, './DragModelsTest/Output/Yu_MeasVsCalc.jpg', 'Resolution', 300);

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
exportgraphics(gcf, './DragModelsTest/Output/Yu_MeasVsCalc_Eqn.jpg', 'Resolution', 300);
