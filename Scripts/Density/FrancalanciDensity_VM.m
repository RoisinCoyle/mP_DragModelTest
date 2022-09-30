%% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Title: Density assumption test: Francalanci
% Date created: 26.06.22
% Date last mostified: 26.06.22
% Purpose: To test the model by Yu satisfies the density and initial
% velocity assumption
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

%% Read in data file
clear
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

%% Francalanci' method 
% <<<<<<<<<<<<<<<<<
% This is not an interative procedure, it just calculates the terminal velocity.

Dref_Frn = zeros(54, 1);
Ddim_Frn = zeros(54, 1);
wdim_Frn = zeros(54, 1);
wvel_Frn = zeros(54, 1);

for i=1:54
    Dref_Frn(i) = size_a(i)*(CSF(i)^0.34)*((size_b(i)/size_a(i))^0.5);
    Ddim_Frn(i) = Dref_Frn(i)*(((g*rho_rel(i))/(vis_kin(i)^2.0))^(1.0/3.0));
	
	E = size_a(i)*((((size_a(i)^2.0)+(size_b(i)^2.0)+(size_c(i)^2.0))/3.0)^-0.5);
	C1 = 18.0*(E^-0.38);
	C2 = 0.3708*(CSF(i)^-0.1602);
	n = 0.4942*(CSF(i)^-0.059);
	
	wdim_Frn(i) = (Ddim_Frn(i)^2)/(C1+(0.75*C2*(Ddim_Frn(i)^3))^n);
	wvel_Frn(i) = wdim_Frn(i)*((rho_rel(i)*g*vis_kin(i))^(1.0/3.0));
end


%% Store output in one array
Results_Frn = zeros(54, 5);

for i=1:54
    Results_Frn(i, 1) = d_equi(i);
    Results_Frn(i, 2) = CSF(i);
    Results_Frn(i, 3) = rho_p(i);
    Results_Frn(i, 4) = wvel_Frn(i);
    Results_Frn(i, 5) = rho_f(i);
end 

% Removing NaN values
n=0;
for i = 1:54
    if (Results_Frn(i, 4) < 0)
        Results_Frn(i, 4) = 0;
    end
end

Table_Frn = array2table(Results_Frn, "VariableNames", ...
    {'ESD', 'CSF', 'rho_p', 'Wt', 'rho_f'});

Table_Frn = [Density_Dataset.ID, Density_Dataset.Shape, Table_Frn];

Table_Frn.Properties.VariableNames(1) = {'ID'};
Table_Frn.Properties.VariableNames(2) = {'Shape'};

writetable(Table_Frn, './DragModelsTest/Output/20220621/Density/FrancalanciOutputDensity.txt', 'Delimiter', ',', 'WriteRowNames', true);
writetable(Table_Frn, './DragModelsTest/Output/20220621/Density/FrancalanciOutputDensity.xls', 'WriteRowNames', true);

%% Plot Frn output
% <<<<<<<<<<<<<<<<<<<
clear
Table_Frn= readtable("./DragModelsTest/Output/20220621/Density/FrancalanciOutputDensity.txt", "Delimiter", ",");

%% Range table

range_m = zeros(9, 1);
rel_difference = zeros(9, 1);
Col_names = ["Shape", "ESD",  "rho_p"];
Var_types = ["cell", "double","double"];
Property_table = table('Size', [9 3], 'VariableTypes', Var_types);
Property_table.Properties.VariableNames = Col_names;

for i=1:6:49
    n= (i+5)/6;
    range_m(n) = (max(Table_Frn.Wt(i:i+5)) - min(Table_Frn.Wt(i:i+5)));
    rel_difference(n) = range_m(n)/(max(Table_Frn.Wt(i:i+5)));
    Property_table(n, 1:3) = [Table_Frn.Shape(i) Table_Frn.ESD(i) Table_Frn.rho_p(i)];
end

range_t = array2table(range_m);
range_t.Properties.VariableNames(1) = {'Wt_range'};

range_output = [Property_table range_t];

writetable(range_output, './DragModelsTest/Output/20220621/Density/FrancalanciRangeWDensity.txt', 'Delimiter', ',');
writetable(range_output, './DragModelsTest/Output/20220621/Density/FrancalanciRangeWDensity.xls', 'WriteRowNames', true);

%% Plot boxplot

label_m = "";
for i=1:54
    if(i<=18)
        label = sprintf('Fragment, %4.1f kg/m3, ESD %4.4fm', Table_Frn.rho_p(i), Table_Frn.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=19 & i<=36)
        label = sprintf('Fibre, %4.1f kg/m3, ESD %4.4fm', Table_Frn.rho_p(i), Table_Frn.ESD(i));
        label_m{i, 1} = label;
    elseif(i>=37 & i<=54)
        label = sprintf('Film, %4.1f kg/m3, ESD %4.4fm', Table_Frn.rho_p(i), Table_Frn.ESD(i));
        label_m{i, 1} = label;
    end
end

label_t = array2table(label_m);
new_table = [Table_Frn label_t];
%%
boxplot(new_table.Wt, new_table.label_m)
ylabel('Terminal Settling Velocity (m/s)')
title(sprintf('Francalanci et al (2021) \n\r %s_f = %5.2f to %5.2f kg/m^3', '\rho', Table_Frn.rho_f(1), Table_Frn.rho_f(6)))
set(gcf, 'WindowState', 'maximized');

exportgraphics(gcf, './DragModelsTest/Output/20220621/Density/Francalanci_Boxplot.jpg', 'Resolution', 300)

%% Calculate 
wvel_list = sort(Table_Frn.Wt(:));
for i=1:6:49
    n= (i+5)/6;
    wvel_range(n) = range(Table_Frn.Wt(i:i+5));
    wvel_iqr(n) = iqr(Table_Frn.Wt(i:i+5));
    wvel_lq(n) = quantile(Table_Frn.Wt(i:i+5), 0.25);
    wvel_uq(n) = quantile(Table_Frn.Wt(i:i+5), 0.75);
    wvel_mid(n) = quantile(Table_Frn.Wt(i:i+5), 0.5);
end

%% A1) wt against Density
% =====================

plot(Table_Frn{1:6, 'rho_f'}, Table_Frn{1:6, 'Wt'}, '-ob', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
hold on
plot(Table_Frn{7:12, 'rho_f'}, Table_Frn{7:12, 'Wt'}, '-sb', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{13:18, 'rho_f'}, Table_Frn{13:18, 'Wt'}, '-^b', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'b')
plot(Table_Frn{19:24, 'rho_f'}, Table_Frn{19:24, 'Wt'}, '-or', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{25:30, 'rho_f'}, Table_Frn{25:30, 'Wt'}, '-sr', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{31:36, 'rho_f'}, Table_Frn{31:36, 'Wt'}, '-^r', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'r')
plot(Table_Frn{37:42, 'rho_f'}, Table_Frn{37:42, 'Wt'}, '-og', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_Frn{43:48, 'rho_f'}, Table_Frn{43:48, 'Wt'}, '-sg', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
plot(Table_Frn{49:54, 'rho_f'}, Table_Frn{49:54, 'Wt'}, '-^g', ...
    'MarkerSize',5,'MarkerEdgeColor','k', 'MarkerFaceColor', 'g')
legend(sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(1), Table_Frn.ESD(1)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(7), Table_Frn.ESD(7)), ...
    sprintf('Fragment, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(13), Table_Frn.ESD(13)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(19), Table_Frn.ESD(19)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(25), Table_Frn.ESD(25)), ...
    sprintf('Fibre, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(31), Table_Frn.ESD(31)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(37), Table_Frn.ESD(37)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(43), Table_Frn.ESD(43)), ...
    sprintf('Film, %4.1f kg/m^{3}, ESD %4.4f m', Table_Frn.rho_p(49), Table_Frn.ESD(49)), ...
    'NumColumns', 3, 'location', 'southoutside')
title('Francalanci et al (2021).')
ylabel('Terminal settling velocity (m/s)')
xlabel('Fluid Density (kg/m^{3})')
   
set(gcf, 'WindowState', 'maximized');
exportgraphics(gcf, './DragModelsTest/Output/20220621/Density/Francalanci_Density.jpg', 'Resolution', 300)