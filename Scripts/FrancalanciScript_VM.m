%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: FrancalanciScript: VM
% Date created: 23.04.22
% Date last mostified: 23.04.22
% Purpose: To test the implementation of the Francalanci drag model on a range of
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

%% Francalanci' method 
% <<<<<<<<<<<<<<<<<
% This is not an interative procedure, it just calculates the terminal velocity.

Dref_Frn = zeros(140, 1);
Ddim_Frn = zeros(140, 1);
wdim_Frn = zeros(140, 1);
wvel_Frn = zeros(140, 1);

for i=1:140
    Dref_Frn(i) = size_a(i)*((CSF(i))^0.34)*((size_b(i)/size_a(i))^0.5);
    Ddim_Frn(i) = Dref_Frn(i)*((g*rho_rel(i))/(vis_kin(i)^2.0))^(1.0/3.0);
	
	E = size_a(i)*((((size_a(i)^2.0)+(size_b(i)^2.0)+(size_c(i)^2.0))/3.0)^-0.5);
	C1 = 18.0*(E^-0.38);
	C2 = 0.3708*(CSF(i)^-0.1602);
	n = 0.4942*(CSF(i)^-0.059);
	
	wdim_Frn(i) = (Ddim_Frn(i)^2)/(C1+(0.75*C2*(Ddim_Frn(i)^3))^n);
	wvel_Frn(i) = ((wdim_Frn(i)*(rho_p(i) - rho_f(i))*g*vis_kin(i))/rho_f(i))^(1.0/3.0);
end

% Store output in one array
Results_Frn = zeros(140, 4);

for i=1:140
    Results_Frn(i, 1) = d_equi(i);
    Results_Frn(i, 2) = CSF(i);
    Results_Frn(i, 3) = wvel_Frn(i);
    Results_Frn(i, 4) = wvel_meas(i);
end 

Table_Frn = array2table(Results_Frn, "VariableNames", ...
    {'ESD', 'CSF', 'Wt','Wt_Meas'});

Table_Frn = [VM_Dataset.Shape Table_Frn];
Table_Frn.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Frn, './DragModelsTest/Output/FrancalanciOutputVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Frn, './DragModelsTest/Output/FrancalanciOutputVM.xls', 'WriteRowNames', true);

%% Plot Francalanci output
% <<<<<<<<<<<<<<<<<<<
Table_Frn= readtable("./DragModelsTest/Output/FrancalanciOutputVM.txt", "Delimiter", ",");

%% A1) wt against ESD
% =====================

plot(Table_Frn.('ESD'), Table_Frn.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Frn.('ESD'), Table_Frn.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Francalanci Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/FrancalanciVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

% Method 1: Shapes Plotted Separately
plot(Table_Frn.('ESD'), Table_Frn.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Frn{1:80, "ESD"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{81:100, "ESD"}, Table_Frn{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "ESD"}, Table_Frn{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Francalanci Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/FrancalanciVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF
% ====================

% Method 1: Plotting all 
plot(Table_Frn.('CSF'), Table_Frn.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Frn.('CSF'), Table_Frn.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Francalanci Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/FrancalanciVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF
% ====================

% Method 1: Shapes Plotted Separately
plot(Table_Frn.('CSF'), Table_Frn.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Frn{1:80, "CSF"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{81:100, "CSF"}, Table_Frn{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "CSF"}, Table_Frn{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Francalanci Model. Using Particle Surface Area')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/FrancalanciVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured
% ============================
Highest(1) = max(Table_Frn.Wt);
Highest(2) = max(Table_Frn.Wt_Meas);
MaxW = max(Highest);
yx=linspace(0, MaxW, 100);

% Method 1: Plot shapes separately
plot(yx, yx)
hold on
plot(Table_Frn{1:80, "Wt_Meas"}, Table_Frn{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{81:100, "Wt_Meas"}, Table_Frn{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{101:140, "Wt_Meas"}, Table_Frn{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Francalanci Model. Using Particle Surface Area')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW*1.1] )
set(gca,'XLim', [0, MaxW*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/Francalanci_MeasVsCalc.jpg', 'Resolution', 300);

%% D) wt against wt measured with fitted lines
% ===============================================

plot(Table_Frn.('Wt_Meas'), Table_Frn.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_Frn.('Wt_Meas'), Table_Frn.('Wt'), 1);
px=[min(Table_Frn.('Wt_Meas')) max(Table_Frn.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(0.75*px(2), py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Frn.("Wt_Meas")\Table_Frn.("Wt");
mx = m*Table_Frn.("Wt_Meas");
plot(Table_Frn.('Wt_Meas'), mx, '-g');
text(0.75*px(2), 0.85*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Francalanci Model. Using Particle Surface Area')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW])
set(gca, 'Xlim', [0, 1.1*MaxW])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/Francalanci_MeasVsCalc_Eqn.jpg', 'Resolution', 300);
