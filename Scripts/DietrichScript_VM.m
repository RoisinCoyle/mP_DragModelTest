%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: DietrichScript: VM
% Date created: 23.04.22
% Date last mostified: 23.04.22
% Purpose: To test the implementation of the Dietrich drag model on a range of
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

%% Dietrich' method 1
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

writetable(Table_Dietrich, './DragModelsTest/Output/DietrichOutputVM.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dietrich, './DragModelsTest/Output/DietrichOutputVM.xls', 'WriteRowNames', true);

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

NewShape = Table_Dietrich.Shape([Particle_Num]);
NewShape_T = array2table(NewShape);
Table_Dietrich_NaN = array2table(Results_NaN_Dietrich, "VariableNames", ...
    {'ESD', 'CSF', 'Wt', 'Wt_Meas'});

Table_Dietrich_NaN = [NewShape_T Table_Dietrich_NaN];
Table_Dietrich_NaN.Properties.VariableNames(1) = {'Shape'};

writetable(Table_Dietrich_NaN, './DragModelsTest/Output/DietrichOutputVM_NaN.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Dietrich_NaN, './DragModelsTest/Output/DietrichOutputVM_NaN.xls', 'WriteRowNames', true);

%% Plot Dietrich output
% <<<<<<<<<<<<<<<<<<<
Table_Dietrich= readtable("./DragModelsTest/Output/DietrichOutputVM.txt", "Delimiter", ",");

%% A1) wt against ESD
% =====================

plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dietrich Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('Particle size (m)')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/DietrichVM_ESDVsW.jpg', 'Resolution', 300)

%% A2) wt against ESD
% =====================

% Method 1: Shapes Plotted Separately
plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
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
exportgraphics(gcf, './DragModelsTest/Output/DietrichVM_ESDVsW_Shapes.jpg', 'Resolution', 300)

%% B1) wt against CSF
% ====================

% Method 1: Plotting all 
plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
legend('Measured Wt', 'Calculated Wt', 'location', 'best')
title('Dietrich Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/DietrichVM_CSFVsW.jpg', 'Resolution', 300);

%% B2) wt against CSF
% ====================

% Method 1: Shapes Plotted Separately
plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dietrich{1:80, "CSF"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich{81:100, "CSF"}, Table_Dietrich{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "CSF"}, Table_Dietrich{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dietrich Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/DietrichVM_CSFVsW_Shapes.jpg', 'Resolution', 300);

%% C) wt against wt measured
% ============================
Highest(1) = max(Table_Dietrich.Wt);
Highest(2) = max(Table_Dietrich.Wt_Meas);
MaxW = max(Highest);
yx=linspace(0, MaxW, 100);

% Method 1: Plot shapes separately
plot(yx, yx)
hold on
plot(Table_Dietrich{1:80, "Wt_Meas"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich{81:100, "Wt_Meas"}, Table_Dietrich{81:100, "Wt"}, 'or',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "Wt_Meas"}, Table_Dietrich{101:140, "Wt"}, 'og',...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
title('Dietrich Model.')
xlabel('Measured Velocity (m/s)')
ylabel('Calculated Velocity (m/s)')
legend('', 'Fragment', 'Fibre', 'Film', 'location', 'best')
set(gca,'YLim', [0, MaxW*1.1] )
set(gca,'XLim', [0, MaxW*1.1] )
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/DietrichVM_MeasVsCalc.jpg', 'Resolution', 300);

%% D) wt against wt measured with fitted lines
% ===============================================

plot(Table_Dietrich.('Wt_Meas'), Table_Dietrich.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_Dietrich.('Wt_Meas'), Table_Dietrich.('Wt'), 1);
px=[min(Table_Dietrich.('Wt_Meas')) max(Table_Dietrich.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 0.9*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Dietrich.("Wt_Meas")\Table_Dietrich.("Wt");
mx = m*Table_Dietrich.("Wt_Meas");
plot(Table_Dietrich.('Wt_Meas'), mx, '-g');
text(px(2), 0.9*max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 7.5, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Dietrich Model.')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW])
set(gca, 'Xlim', [0, 1.1*MaxW])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/DietrichVM_MeasVsCalc_Eqn.jpg', 'Resolution', 300);

%% Plot output after removing NaN Values
Table_Dietrich_New = readtable('./DragModelsTest/Output/DietrichOutputVM_NaN.txt', 'Delimiter', ',');

Highest(1) = max(Table_Dietrich_New.Wt);
Highest(2) = max(Table_Dietrich_New.Wt_Meas);
MaxW = max(Highest);
yx=linspace(0, MaxW, 100);

%%
plot(Table_Dietrich_New.('Wt_Meas'), Table_Dietrich_New.('Wt'), 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(yx, yx, '-k')
p=polyfit(Table_Dietrich_New.('Wt_Meas'), Table_Dietrich_New.('Wt'), 1);
px=[min(Table_Dietrich_New.('Wt_Meas')) max(Table_Dietrich_New.('Wt_Meas'))];
py=polyval(p, px);
plot(px, py, '-b')
text(px(2), 1.05*py(2), (sprintf('y = %.4fx %+.4f', p(1), p(2))), ...
    'Color', 'b', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
m=Table_Dietrich_New.("Wt_Meas")\Table_Dietrich_New.("Wt");
mx = m*Table_Dietrich_New.("Wt_Meas");
plot(Table_Dietrich_New.('Wt_Meas'), mx, '-g');
text(px(2), max(mx), (sprintf('y = %.4fx', m)), ...
    'Color', 'g', 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left');
title('Dietrich Model.')
xlabel('Measured Wt (m/s)')
ylabel('Calculated Wt (m/s)')
legend('', 'y=x', 'Linear fit', 'Linear fit forced', 'location', 'best')
set(gca, 'Ylim', [0, 1.1*MaxW])
set(gca, 'Xlim', [0, 1.1*MaxW])
hold off

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/DietrichVM_MeasVsCalc_Eqn_NaN.jpg', 'Resolution', 300);

%%
subplot(1, 2, 1)
plot(Table_Dietrich.('CSF'), Table_Dietrich.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
hold on
plot(Table_Dietrich{1:80, "CSF"}, Table_Dietrich{1:80, "Wt"}, 'ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Dietrich{81:100, "CSF"}, Table_Dietrich{81:100, "Wt"}, 'or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Dietrich{101:140, "CSF"}, Table_Dietrich{101:140, "Wt"}, 'og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend('Measured Wt', 'Calculated Wt, Fragment', 'Calculated Wt, Fibre', ...
       'Calculated Wt, Film', 'NumColumns', 2, 'location', 'southoutside')
title('Dietrich Model.')
ylabel('Terminal settling velocity (m/s)')
xlabel('CSF')
hold off

subplot(1, 2, 2)
plot(Table_Dietrich.('ESD'), Table_Dietrich.('Wt_Meas'), 'ok', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'k')
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
exportgraphics(gcf, './DragModelsTest/Output/DietrichVM_ESDV_CSF_NaN.jpg', 'Resolution', 300)