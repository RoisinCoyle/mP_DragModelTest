%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: YuScript: VM
% Date created: 25.04.22
% Date last mostified: 02.03.23
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
Results_Yu = zeros(54, 5);

for i=1:54
    Results_Yu(i, 1) = d_equi(i);
    Results_Yu(i, 2) = CSF(i);
    Results_Yu(i, 3) = wvel_Yu(i);
    Results_Yu(i, 4) = rho_f(i);
    Results_Yu(i, 5) = rho_p(i);
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
    {'ESD', 'CSF', 'Wt', 'rho_f', 'rho_p'});

Table_Yu = [Density_Dataset.ID, Density_Dataset.Shape, Table_Yu];

Table_Yu.Properties.VariableNames(1) = {'ID'};
Table_Yu.Properties.VariableNames(2) = {'Shape'};

writetable(Table_Yu, './DragModelsTest/Output/20220621/Density/YuOutputDensity.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Yu, './DragModelsTest/Output/20220621/Density/YuOutputDensity.xls', 'WriteRowNames', true);

%% Read in data
% <<<<<<<<<<<<<<<<<<<
clear
Table_Yu= readtable("./DragModelsTest/Output/20220621/Density/YuOutputDensity.txt", "Delimiter", ",");

label_m = "";
for i=1:54
    if(i<=18)
        label = sprintf('Fragment, %4.1f kg/m3, ESD %4.4fm', Table_Yu.rho_p(i), Table_Yu.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=19 & i<=36)
        label = sprintf('Fibre, %4.1f kg/m3, ESD %4.4fm', Table_Yu.rho_p(i), Table_Yu.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=37 & i<=54)
        label = sprintf('Film, %4.1f kg/m3, ESD %4.4fm', Table_Yu.rho_p(i), Table_Yu.ESD(i));
        label_m{i, 1} = label;
    end
end

label_t = array2table(label_m);
new_table = [Table_Yu label_t];
%% Plot boxplot
colors = {[0.6980 1 0.4118] [0.6980 1 0.4118] [0.6980 1 0.4118] ...
          [1 0.6000 0.6000] [1 0.6000 0.6000] [1 0.6000 0.6000] ...
          [0.4000 0.6980 1] [0.4000 0.6980 1] [0.4000 0.6980 1] };

fig = figure
hold on
boxplot(new_table.Wt, new_table.label_m, 'position', (1:6:54), 'widths', 5, 'boxstyle', 'outline', 'Colors', 'k')
ylim([0 0.035])
boxes = fig.Children.Children(1,1).Children(19:27)
for j = 1:length(boxes) % draw a colored patch behind each bar
        patch(boxes(j).XData, boxes(j).YData, colors{j},'FaceAlpha',.5,'EdgeAlpha',0.3);
end
ylabel('Modelled Terminal Settling Velocity (m/s)')
title(sprintf('Boxplots showing the range of modelled terminal settling velocity attained when fluid density varies from %s_f = %5.2f to %5.2f kg/m^3.', '\rho', Table_Yu.rho_f(1), Table_Yu.rho_f(6)))
subtitle('Model applied: Yu et al (2022).')

set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Density/Yu_Boxplot.jpg', 'Resolution', 1200)

%% Range table

range_m = zeros(9, 1);
rel_difference = zeros(9, 1);
Col_names = ["Shape", "ESD",  "rho_p"];
Var_types = ["cell", "double","double"];
Property_table = table('Size', [9 3], 'VariableTypes', Var_types);
Property_table.Properties.VariableNames = Col_names;

for i=1:6:49
    n= (i+5)/6;
    range_m(n) = (max(Table_Yu.Wt(i:i+5)) - min(Table_Yu.Wt(i:i+5)));
    rel_difference(n) = range_m(n)/(max(Table_Yu.Wt(i:i+5)));
    Property_table(n, 1:3) = [Table_Yu.Shape(i) Table_Yu.ESD(i) Table_Yu.rho_p(i)];
end

range_t = array2table(range_m);
range_t.Properties.VariableNames(1) = {'Wt_range'};

range_output = [Property_table range_t];

writetable(range_output, './DragModelsTest/Output/20220621/Density/YuRangeWDensity.txt', 'Delimiter', ',');
writetable(range_output, './DragModelsTest/Output/20220621/Density/YuRangeWDensity.xls', 'WriteRowNames', true);
%% A1) wt against Density
% =====================

plot(Table_Yu{1:6, 'rho_f'}, Table_Yu{1:6, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Yu{7:12, 'rho_f'}, Table_Yu{7:12, 'Wt'}, '-sb', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{13:18, 'rho_f'}, Table_Yu{13:18, 'Wt'}, '-^b', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Yu{19:24, 'rho_f'}, Table_Yu{19:24, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{25:30, 'rho_f'}, Table_Yu{25:30, 'Wt'}, '-sr', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{31:36, 'rho_f'}, Table_Yu{31:36, 'Wt'}, '-^r', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Yu{37:42, 'rho_f'}, Table_Yu{37:42, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_Yu{43:48, 'rho_f'}, Table_Yu{43:48, 'Wt'}, '-sg', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_Yu{49:54, 'rho_f'}, Table_Yu{49:54, 'Wt'}, '-^g', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend(sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(1), Table_Yu.ESD(1)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(7), Table_Yu.ESD(7)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(13), Table_Yu.ESD(13)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(19), Table_Yu.ESD(19)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(25), Table_Yu.ESD(25)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(31), Table_Yu.ESD(31)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(37), Table_Yu.ESD(37)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(43), Table_Yu.ESD(43)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Yu.rho_p(49), Table_Yu.ESD(49)), ...
    'NumColumns', 3, 'location', 'southoutside')
title(sprintf("The impact of fluid density (%s_f) on modelled terminal settling velocity of six mP particles selected randomly from Van Melkebeke et al (2020)'s dataset.", '\rho'))
subtitle('Model applied: Yu et al (2022).')
ylabel('Modelled Terminal settling velocity (m/s)')
xlabel('Fluid Density (kg/m^{3})')

   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20230301/Density/Yu_Density.jpg', 'Resolution', 1200)
